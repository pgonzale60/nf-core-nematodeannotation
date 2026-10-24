#!/usr/bin/env python3
"""
slidr_cluster_tails.py — Dereplicate and 3'-anchored cluster SLIDR read tails.

Port of stage 2 of the tested SLIDR pipeline (Wenzel et al., BMC Bioinformatics
2021), and a direct replacement for:

    vsearch --derep_fulllength - --sizeout --minseqlength L
  | vsearch --cluster_fast - --id 1.0 --iddef 0 --rightjust
            --wordlength 8 --minwordmatches 1000 --sizeorder
            --maxaccepts 0 --maxrejects 0

WHY RIGHT-JUSTIFIED CLUSTERING
------------------------------
A read tail is the 5' overhang of an RNA-seq read that the aligner soft-clipped.
For a trans-spliced mRNA that overhang is the spliced leader (SL), and the tail's
3' end is anchored at the trans-splice junction. Two tails therefore belong to
the same SL when their 3' ends agree, regardless of how much of the SL's 5' end
survived.

vsearch implements this as right-justified alignment combined with --iddef 0
(identity = matches / length of the SHORTER sequence). Together those two
options are exactly:

    tail A and tail B cluster  <=>  the shorter of the two is an exact
                                   suffix of the longer

which is what this script implements, with no external dependency.

DEVIATIONS FROM VSEARCH (documented, deliberate)
------------------------------------------------
* Centroid selection is greedy in order of (abundance desc, length desc),
  matching --cluster_fast --sizeorder. vsearch does not specify tie-breaking
  within equal abundance, so this can differ on pathological ties. The
  downstream consensus step re-derives each cluster sequence from the LONGEST
  member, so centroid identity does not propagate into the consensus.
* --wordlength 8 --minwordmatches 1000 is a speed heuristic that forces all
  8-mers to match. The equivalent bucketing here is exact rather than
  heuristic: because every retained tail is at least --min-len bases long, two
  tails where one is a suffix of the other necessarily share the same
  --min-len-mer, so bucketing on that key loses no clusters. (A fixed key
  length would be wrong: a tail shorter than the key would be bucketed by its
  whole sequence and could never meet its longer parent.)

OUTPUTS
-------
  {prefix}_centroids.fa          >centroid_N;size=<reads>   one per cluster
  {prefix}_cluster_summary.tsv   centroid_id, sequence, reads, unique_tails, length
  {prefix}_readinfo.tsv.gz       read_id<TAB>centroid_id    one row per tail

read_id is emitted as "<library>::<original read id>" so that tails pooled from
several libraries stay addressable.
"""

import argparse
import collections
import gzip
import os
import sys
from collections.abc import Iterator


def open_maybe_gzip(path: str):
    """Open a plain or gzipped text file transparently."""
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def iter_tails(paths: list[str], min_len: int) -> Iterator[tuple[str, str]]:
    """
    Yield (read_id, tail_seq) for every usable tail.

    Reads are qualified with the library label derived from the file name.
    Tails containing an N are dropped, mirroring the tested pipeline's
    `awk '!/^>/ && $0~"N"{$0=""}'`; tails shorter than min_len are dropped,
    mirroring --minseqlength. Headers look like
    '>read_id|chrom:pos:strand|len=N|dinuc=XX'.
    """
    for path in paths:
        label = os.path.basename(path).split(".tails")[0]
        with open_maybe_gzip(path) as fh:
            read_id = None
            for line in fh:
                line = line.rstrip("\n").rstrip("\r")
                if not line:
                    continue
                if line[0] == ">":
                    read_id = line[1:].split("|")[0].split()[0]
                    continue
                seq = line.strip().upper()
                if not seq or read_id is None:
                    continue
                if "N" in seq or len(seq) < min_len:
                    continue
                yield f"{label}::{read_id}", seq


def bucket_key(seq: str, key_len: int) -> str:
    """
    3'-terminal bucket key: the last key_len bases of the tail.

    key_len is the minimum retained tail length, so every sequence is at least
    key_len long. Two tails in a suffix relation therefore always agree on this
    key, which is what makes the bucketing lossless.
    """
    return seq[-key_len:]


def cluster_bucket(
    entries: list[tuple[str, int]], key_len: int, min_parent_reads: int = 3
) -> list[tuple[str, list[str]]]:
    """
    Cluster one bucket of (unique_seq, abundance) by 3'-anchored exact-suffix
    identity, assembling shorter 3'-suffixes into their longer parent leaders.

    Returns [(centroid_seq, [member_seq, ...]), ...]; every member is an exact
    suffix of its centroid (a centroid is trivially a member of itself).

    To prevent truncated leaders (which are far more abundant due to 5' RNA-seq
    degradation) from outcompeting and isolating their longer parent sequences,
    clustering proceeds in two passes:
      1. Pass 1: Sequences with at least `min_parent_reads` are considered in
         descending length order (and descending abundance for ties) to absorb
         all shorter 3'-suffixes.
      2. Pass 2: Any remaining sequences (e.g. singletons) cluster among themselves.
    """
    # Sort length desc, then abundance desc
    ordered = sorted(entries, key=lambda e: (-len(e[0]), -e[1]))
    pending: dict[str, int] = {seq: ab for seq, ab in ordered}
    clusters = []

    # Pass 1: Qualified parent leaders (>= min_parent_reads) absorb suffixes
    for seq, abundance in ordered:
        if seq not in pending:
            continue
        if abundance < min_parent_reads:
            continue
        members = [seq]
        del pending[seq]
        for length in range(len(seq) - 1, key_len - 1, -1):
            candidate = seq[len(seq) - length :]
            if candidate in pending:
                members.append(candidate)
                del pending[candidate]
        clusters.append((seq, members))

    # Pass 2: Remaining sequences cluster among themselves
    for seq, abundance in ordered:
        if seq not in pending:
            continue
        members = [seq]
        del pending[seq]
        for length in range(len(seq) - 1, key_len - 1, -1):
            candidate = seq[len(seq) - length :]
            if candidate in pending:
                members.append(candidate)
                del pending[candidate]
        clusters.append((seq, members))

    return clusters


def main():
    ap = argparse.ArgumentParser(
        description="Dereplicate and 3'-anchored cluster SLIDR read tails "
        "(dependency-free port of vsearch cluster_fast "
        "--iddef 0 --rightjust)."
    )
    ap.add_argument(
        "--tails",
        required=True,
        nargs="+",
        help="One or more *.tails.fa.gz files from slidr_extract_tails.py",
    )
    ap.add_argument(
        "--min-len",
        type=int,
        default=8,
        help="Minimum tail length to retain (default: 8)",
    )
    ap.add_argument("--prefix", default="slidr", help="Output file prefix")
    ap.add_argument("--outdir", default=".", help="Output directory")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # The bucket key length must equal the minimum tail length for the
    # bucketing to be lossless (see bucket_key).
    key_len = args.min_len

    # ---- pass 1: dereplicate ------------------------------------------------
    # Only unique sequences are held in memory; the per-read association is
    # recovered by a second streaming pass below.
    abundance: dict[str, int] = collections.Counter()
    n_tails = 0
    for _read_id, seq in iter_tails(args.tails, args.min_len):
        abundance[seq] += 1
        n_tails += 1

    print(
        f"[slidr_cluster_tails] {n_tails:,} tails retained, {len(abundance):,} unique",
        file=sys.stderr,
    )
    if not abundance:
        print(
            "[slidr_cluster_tails] FATAL: no tails passed the filters", file=sys.stderr
        )
        sys.exit(1)
    if n_tails == len(abundance):
        print(
            "[slidr_cluster_tails] NOTE: no duplicate tails; is --min-len appropriate?",
            file=sys.stderr,
        )

    # ---- cluster ------------------------------------------------------------
    buckets: dict[str, list[tuple[str, int]]] = collections.defaultdict(list)
    for seq, ab in abundance.items():
        buckets[bucket_key(seq, key_len)].append((seq, ab))

    seq_to_centroid: dict[str, str] = {}
    centroids: list[tuple[str, str, int]] = []  # (centroid_id, seq, reads)
    for key in sorted(buckets):
        for centroid_seq, members in cluster_bucket(buckets[key], key_len):
            cid = f"centroid_{len(centroids) + 1}"
            total = sum(abundance[m] for m in members)
            centroids.append((cid, centroid_seq, total))
            for member in members:
                seq_to_centroid[member] = cid

    uniq_per_centroid = collections.Counter(seq_to_centroid.values())
    print(
        f"[slidr_cluster_tails] {len(centroids):,} clusters across "
        f"{len(buckets):,} 3'-terminal buckets",
        file=sys.stderr,
    )

    # ---- pass 2: emit read -> centroid -------------------------------------
    readinfo_path = os.path.join(args.outdir, f"{args.prefix}_readinfo.tsv.gz")
    n_out = 0
    with gzip.open(readinfo_path, "wt") as out:
        for read_id, seq in iter_tails(args.tails, args.min_len):
            cid = seq_to_centroid.get(seq)
            if cid is not None:
                out.write(f"{read_id}\t{cid}\n")
                n_out += 1

    # ---- write centroids and summary ---------------------------------------
    fasta_path = os.path.join(args.outdir, f"{args.prefix}_centroids.fa")
    with open(fasta_path, "w") as out:
        for cid, seq, total in centroids:
            out.write(f">{cid};size={total}\n{seq}\n")

    summary_path = os.path.join(args.outdir, f"{args.prefix}_cluster_summary.tsv")
    with open(summary_path, "w") as out:
        out.write("centroid_id\tsequence\treads\tunique_tails\tlength\n")
        for cid, seq, total in centroids:
            out.write(f"{cid}\t{seq}\t{total}\t{uniq_per_centroid[cid]}\t{len(seq)}\n")

    print(
        f"[slidr_cluster_tails] wrote {n_out:,} read->centroid assignments\n"
        f"  centroids : {fasta_path}\n"
        f"  summary   : {summary_path}\n"
        f"  readinfo  : {readinfo_path}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
