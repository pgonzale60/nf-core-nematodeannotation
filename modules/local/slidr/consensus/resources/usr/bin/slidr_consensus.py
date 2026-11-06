#!/usr/bin/env python3
"""
slidr_consensus.py — Resolve consensus spliced leaders and trans-splice acceptor sites.

This replaces the earlier, simplified port with the logic of the tested SLIDR
consensus stage (slidr_consensus.R, Wenzel et al., BMC Bioinformatics 2021). The
simplified version decided leader identity by COMPARING TAIL SEQUENCES TO A
HARD-CODED TABLE of known leaders, which cannot work for an organism whose
leaders are not in that table -- and Auanema rhodense is such an organism: the
canonical SL2 leader occurs zero times in its assembly and zero times at read 5'
ends, while SL1 is abundant.

Here a tail is attributed to a leader by JOINING three tables produced upstream:

    read --(SLIDR_CLUSTERTAILS)--> tail cluster centroid
    centroid --(SLIDR_SLRNASITES)--> leader RNA gene locus (3' end aligns to the
                                genome, splice donor within max-overlap bases,
                                Sm-binding motif within max-rna-len bases)

Leaders are then NAMED DE NOVO, SL1, SL2, ... in order of the number of
trans-splice sites they support. Nothing about the answer is assumed in advance.

Filters retained from the tested implementation
-----------------------------------------------
1. Acceptor motif. A site is only kept when the genomic dinucleotide immediately
   upstream of the alignment start matches the acceptor pattern (AG by default).
   SLIDR_EXTRACTTAILS has already evaluated this per read (is_canonical_ag).

2. Leader-gene distance. A trans-spliced gene's acceptor site must be further than
   --max-leader-gene-distance bp from the leader RNA gene that served it. Without
   this the leader tandem array's own transcripts (trivially spliced to the
   leader) swamp the real sites. Sites on a different chromosome or strand from
   their leader are kept, as distance is undefined there.

3. Cluster consensus by majority length. Within a cluster the consensus is the
   longest member, and a member counts as full length only if it is within the
   median length deficit of that longest member. This is how the tested
   implementation separates a genuine leader from 5'-truncated fragments of it --
   which matters here, because the leaders in these libraries arrive as ~11 nt 3'
   fragments rather than intact 22-mers.

Only leaders with at least --min-support reads AND at least 2 independent sites
are reported, mirroring the tested "at least 2 reads and spliced to at least 2
sites" rule.

Outputs (names are a hard contract with Agent 3; do not rename):
    {prefix}_sl_acceptor_sites.bed    BED6
    {prefix}_sl_acceptor_sites.gff3   GFF3
    {prefix}_sl_summary.tsv           per-leader summary
    {prefix}_sl_sequences.fa          consensus leader sequences
"""

import argparse
import collections
import gzip
import os
import re
import statistics
import sys
from typing import Dict, List, Optional, Tuple


def open_maybe_gzip(path: str):
    """Open a plain or gzipped text file transparently."""
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def library_label(path: str) -> str:
    """
    Library label for a per-library file, used to key reads.

    SLIDR_CLUSTERTAILS derives its label from '<library>.tails.fa.gz'; the site
    tables here are '<library>.acceptor_sites.tsv.gz'. Both reduce to the same
    '<library>' prefix, which is what lets the join work.
    """
    base = os.path.basename(path)
    for marker in (".acceptor_sites", ".tails"):
        if marker in base:
            return base.split(marker)[0]
    return base.split(".")[0]


def load_readinfo(path: str) -> Dict[str, str]:
    """
    read_id -> centroid_id, from SLIDR_CLUSTERTAILS.

    Reads absent from this map are those whose tail never reached clustering (too
    short, or containing an N); they are dropped from the analysis.
    """
    mapping = {}
    with open_maybe_gzip(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                mapping[parts[0]] = parts[1]
    return mapping


def load_leader_loci(path: str) -> Dict[str, List[dict]]:
    """
    centroid_id -> [leader RNA locus, ...], from SLIDR_SLRNASITES.

    Columns: centroid_id, chrom, sl_start, sl_end, strand, sl_seq, donor_motif,
             sm_motif, sl_rna_seq
    """
    loci: Dict[str, List[dict]] = collections.defaultdict(list)
    with open(path, "rt") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[0] == "centroid_id":
                continue
            try:
                loci[parts[0]].append({
                    "chrom": parts[1],
                    "start": int(parts[2]),
                    "end": int(parts[3]),
                    "strand": parts[4],
                    "sl_seq": parts[5],
                    "sl_rna_seq": parts[8],
                })
            except ValueError:
                continue
    return loci


def passes_leader_distance(pos: int, chrom: str, strand: str, locus: dict,
                           max_dist: int) -> bool:
    """
    True when the acceptor site is far enough from the leader RNA gene to be a
    genuine trans-splice site of some other gene.

    A leader gene is transcribed as part of the leader array, so reads from that
    array are spliced to the leader trivially, without any trans-splicing. Those
    sites sit immediately next to the leader locus and must be discarded, or they
    dominate the output. Distance is undefined across chromosomes/strands, where
    the testing implementation also keeps the site.
    """
    if chrom != locus["chrom"] or strand != locus["strand"]:
        return True
    if pos < locus["start"]:
        d = locus["start"] - pos
    elif pos > locus["end"]:
        d = pos - locus["end"]
    else:
        d = 0
    return d > max_dist


def load_rrna_intervals(path: Optional[str], flank: int = 500) -> Dict[str, List[Tuple[int, int]]]:
    """
    Load rRNA/rDNA feature intervals from GFF3 (Barrnap or Infernal cmscan output).
    If a comprehensive ncRNA GFF3 is supplied, filters specifically for rRNA features
    (e.g., Name/Accession matching rRNA or canonical Rfam rRNA models RF00001, RF00002, etc.).
    Expands each interval by +/- flank bp to catch reads mapping to intergenic
    spacers or unannotated elements within the rDNA repeat arrays.
    """
    if not path or not os.path.exists(path):
        return {}
    intervals = collections.defaultdict(list)
    rrna_patterns = re.compile(r"rRNA|ribosomal|RF00001|RF00002|RF01959|RF01960|RF02540|RF02541|RF02542|RF02543", re.IGNORECASE)
    with open_maybe_gzip(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            # If line has attributes, ensure it represents an rRNA feature
            if len(parts) >= 9 and not rrna_patterns.search(parts[8]) and not rrna_patterns.search(parts[2]):
                continue
            chrom = parts[0]
            try:
                start = max(0, int(parts[3]) - 1 - flank)
                end = int(parts[4]) + flank
                intervals[chrom].append((start, end))
            except ValueError:
                continue
    # Merge overlapping intervals per chromosome for fast binary/linear search
    merged = {}
    for chrom, ivs in intervals.items():
        ivs.sort()
        m = []
        for s, e in ivs:
            if not m or s > m[-1][1]:
                m.append([s, e])
            else:
                m[-1][1] = max(m[-1][1], e)
        merged[chrom] = m
    return merged


def is_in_rrna(chrom: str, pos_0based: int, rrna_intervals: Dict[str, List[List[int]]]) -> bool:
    """Check if a coordinate falls inside any expanded rRNA/rDNA interval."""
    if chrom not in rrna_intervals:
        return False
    for s, e in rrna_intervals[chrom]:
        if s <= pos_0based <= e:
            return True
        if s > pos_0based:
            break
    return False


def load_sites(paths: List[str], readinfo: Dict[str, str],
               leaders: Dict[str, List[dict]], acceptor: str, max_dist: int,
               rrna_intervals: Optional[Dict[str, List[List[int]]]] = None):
    """
    Stream every per-library site table, keeping only reads that can be
    attributed to a discovered leader, excluding rRNA/rDNA mapping reads.

    Returns (sites, stats) where each site is
    {read, centroid, leader_key, chrom, pos, strand, tail}.
    """
    sites = []
    stats = collections.Counter()
    if rrna_intervals is None:
        rrna_intervals = {}

    for path in paths:
        lib = library_label(path)
        with open_maybe_gzip(path) as fh:
            header = fh.readline().rstrip("\n").split("\t")
            idx = {name: i for i, name in enumerate(header)}
            required = ("read_id", "chrom", "junction_pos_1based", "strand",
                        "tail_seq", "genomic_dinuc", "is_canonical_ag")
            missing = [c for c in required if c not in idx]
            if missing:
                print(f"[slidr_consensus] WARNING: {path} is missing {missing}; skipped",
                      file=sys.stderr)
                continue
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) <= max(idx.values()):
                    continue
                stats["rows"] += 1

                try:
                    pos = int(parts[idx["junction_pos_1based"]])
                except ValueError:
                    continue
                chrom = parts[idx["chrom"]]
                strand = parts[idx["strand"]]

                # Exclude reads falling in annotated rRNA/rDNA arrays or their flanks
                if is_in_rrna(chrom, pos - 1, rrna_intervals):
                    stats["dropped_rrna"] += 1
                    continue
                if parts[idx["is_canonical_ag"]].lower() != "true":
                    stats["dropped_non_canonical"] += 1
                    continue
                # Re-apply the acceptor pattern to the reported dinucleotide, as
                # the tested implementation does, so a non-AG acceptor is honoured
                # if the run is configured for one.
                dinuc = parts[idx["genomic_dinuc"]].upper()
                if not re.fullmatch(acceptor, dinuc):
                    stats["dropped_acceptor_mismatch"] += 1
                    continue
                centroid = readinfo.get(f"{lib}::{parts[idx['read_id']]}")
                if centroid is None:
                    stats["dropped_unclustered"] += 1
                    continue
                loci = leaders.get(centroid)
                if not loci:
                    stats["dropped_no_leader"] += 1
                    continue
                try:
                    pos = int(parts[idx["junction_pos_1based"]])
                except ValueError:
                    continue
                chrom = parts[idx["chrom"]]
                strand = parts[idx["strand"]]
                relevant = [l for l in loci
                            if passes_leader_distance(pos, chrom, strand, l, max_dist)]
                if not relevant:
                    stats["dropped_leader_array"] += 1
                    continue
                stats["kept"] += 1
                sites.append({
                    "read": f"{lib}::{parts[idx['read_id']]}",
                    "centroid": centroid,
                    # A centroid may be supported by several leader loci; the
                    # locus listed first is the best-scoring BLASTN hit.
                    "leader_key": centroid,
                    "chrom": chrom,
                    "pos": pos,
                    "strand": strand,
                    "tail": parts[idx["tail_seq"]].upper(),
                })
    return sites, stats


def cluster_consensus(tails: List[str]) -> Tuple[str, int]:
    """
    Majority-length consensus, ported from the tested implementation.

    The consensus is the longest member of the cluster with its first `cutoff`
    bases lower-cased, where cutoff is the median amount by which members fall
    short of the longest member. Lower case therefore marks the 5' stretch that
    the majority of reads do not reach -- which is exactly the part of a leader
    that cannot be trusted in a library where leaders arrive truncated.
    """
    lengths = [len(t) for t in tails]
    longest = max(lengths)
    deficits = sorted(set(longest - l for l in lengths))
    cutoff = int(statistics.median(deficits)) if deficits else 0
    consensus = next(t for t in tails if len(t) == longest)
    masked = consensus[:cutoff].lower() + consensus[cutoff:]
    return masked, cutoff


def is_adapter_contaminated(seq: str) -> bool:
    """Check if sequence contains common Illumina/sequencing adapter artifacts."""
    seq_up = seq.upper()
    adapter_patterns = [
        r"AGATCGGAAG",
        r"GATCGGAAG",
        r"TCGGAAGAG",
        r"CTGTCTCTTATA",
        r"TGGAATTCTCGG",
    ]
    for pat in adapter_patterns:
        if re.search(pat, seq_up):
            return True
    return False


def classify_sl_family(seq: str) -> str:
    """
    Classify spliced leader sequence into biological families based on 3' conserved core motif.
    Mirrors validated empirical benchmarks across nematode taxa:
      - SL1_canonical: Nematode SL1 (e.g. GGTTTAATTACCCAAGTTTG / AATTACCCAAGT...)
      - SL2_canonical_like: Classical SL2 (e.g. GGTTTTAACCCAGTTACTCA / GGTTTTTACCCAGTATCTCA)
      - SL_operon_variant: Clade-specific operon leaders (e.g. ...[AT]AGAGTTTTTG$ / ...GATTTTTTG$)
      - PolyA_PolyT_artifact: Homopolymer soft-clip artifacts
      - Other_minor: Heterogeneous / novel leader candidates
    """
    seq_up = seq.upper()
    if "AATTACCCAAGTTTG" in seq_up or "AATTACCCAAGTTAG" in seq_up or "TTACCCAAGTTTG" in seq_up:
        return "SL1_canonical"
    elif "TTACCCAGTATCTCA" in seq_up or "AACCCAGT" in seq_up or "TTACCCAGT" in seq_up or "AACCCAAACTT" in seq_up:
        return "SL2_canonical_like"
    elif re.search(r"[AT]AGAGTTTTTG$", seq_up) or re.search(r"GATTTTTTG$", seq_up):
        return "SL_operon_variant"
    elif re.search(r"T{8,}[GATC]?$", seq_up) or re.search(r"A{8,}[GATC]?$", seq_up):
        return "PolyA_PolyT_artifact"
    else:
        return "Other_minor"


def main():
    ap = argparse.ArgumentParser(
        description="Resolve consensus spliced leaders and trans-splice acceptor "
                    "sites by joining clustered tails to de-novo leader loci.")
    ap.add_argument("--sites", required=True, nargs="+",
                    help="Per-library *.acceptor_sites.tsv.gz from SLIDR_EXTRACTTAILS")
    ap.add_argument("--readinfo", required=True,
                    help="*_readinfo.tsv.gz (read -> centroid) from SLIDR_CLUSTERTAILS")
    ap.add_argument("--slrnas", required=True,
                    help="*.sl_rna_candidates.tsv from SLIDR_SLRNASITES")
    ap.add_argument("--acceptor", default="AG",
                    help="Acceptor site pattern, regex (default: AG)")
    ap.add_argument("--min-support", type=int, default=5,
                    help="Minimum supporting reads per reported acceptor site (default: 5)")
    ap.add_argument("--consolidation-len", type=int, default=15,
                    help="Length of 3'-terminal k-mer to consolidate fragmented centroids sharing a leader (default: 15, 0 to disable)")
    ap.add_argument("--max-leader-gene-distance", type=int, default=1000,
                    help="Discard sites within this many bp of their own leader RNA "
                         "gene, i.e. the leader array's own transcripts (default: 1000)")
    ap.add_argument("--rrna-gff", default=None,
                    help="GFF file of annotated rRNA/rDNA genes (e.g. from Barrnap) to exclude")
    ap.add_argument("--rrna-flank", type=int, default=500,
                    help="Flanking distance (bp) around rRNA features to exclude (default: 500)")
    ap.add_argument("--rfam-gff", default=None,
                    help="Optional ncRNA GFF3 (e.g. Infernal cmscan) to cross-reference SL loci")
    ap.add_argument("--prefix", default="aps4", help="Output file prefix")
    ap.add_argument("--outdir", default=".", help="Output directory")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    readinfo = load_readinfo(args.readinfo)
    leaders = load_leader_loci(args.slrnas)
    rrna_intervals = load_rrna_intervals(args.rrna_gff, args.rrna_flank)
    rfam_loci = load_rfam_sl_loci(args.rfam_gff)
    if rrna_intervals:
        n_rrna_regions = sum(len(ivs) for ivs in rrna_intervals.values())
        print(f"[slidr_consensus] Loaded {n_rrna_regions} expanded rRNA/rDNA exclusion zones "
              f"from {args.rrna_gff} (+/- {args.rrna_flank} bp)", file=sys.stderr)
    if rfam_loci:
        print(f"[slidr_consensus] Loaded {len(rfam_loci)} Rfam SL loci for orthogonal validation "
              f"from {args.rfam_gff}", file=sys.stderr)

    print(f"[slidr_consensus] {len(readinfo):,} clustered reads, "
          f"{len(leaders):,} centroids with a leader RNA locus", file=sys.stderr)

    sites, stats = load_sites(args.sites, readinfo, leaders, args.acceptor,
                              args.max_leader_gene_distance, rrna_intervals)
    print("[slidr_consensus] site filtering:", file=sys.stderr)
    for key in ("rows", "dropped_rrna", "dropped_non_canonical", "dropped_acceptor_mismatch",
                "dropped_unclustered", "dropped_no_leader", "dropped_leader_array",
                "kept"):
        print(f"    {key:<28} {stats[key]:>12,}", file=sys.stderr)

    if not sites:
        print("[slidr_consensus] FATAL: no sites survived filtering. Check that "
              "SLIDR_ALIGNREADS actually soft-clipped leader sequence: if no "
              "centroid has a leader RNA locus, the aligner was not configured to "
              "clip (see modules/sanger-tol/slidr/alignreads/main.nf).",
              file=sys.stderr)
        sys.exit(1)

    # ---- locus-level pool (no per-centroid threshold) ------------------------
    # A single leader's reads are split by clustering across many centroids, so
    # the per-(centroid, site) --min-support below under-counts real sites. This
    # table pools every leader-attributed read at the genomic site instead, and
    # leaves thresholding to the consumer (SLIDR_TSSEXPORTER).
    pool = collections.defaultdict(lambda: {"reads": set(), "libs": set(),
                                            "centroids": set(),
                                            "tails": collections.Counter()})
    for s in sites:
        p = pool[(s["chrom"], s["pos"], s["strand"])]
        p["reads"].add(s["read"])
        p["libs"].add(s["read"].split("::", 1)[0])
        p["centroids"].add(s["centroid"])
        p["tails"][s["tail"]] += 1
    pool_path = os.path.join(args.outdir, f"{args.prefix}_sl_site_pool.tsv")
    with open(pool_path, "w") as out:
        out.write("chrom\tpos\tstrand\treads\tn_libraries\tn_centroids\t"
                  "top_tail\ttop_tail_reads\n")
        for (chrom, pos, strand), p in sorted(pool.items()):
            top_tail, top_n = p["tails"].most_common(1)[0]
            out.write(f"{chrom}\t{pos}\t{strand}\t{len(p['reads'])}\t{len(p['libs'])}\t"
                      f"{len(p['centroids'])}\t{top_tail}\t{top_n}\n")
    print(f"[slidr_consensus] {len(pool):,} pooled acceptor loci written to "
          f"{pool_path}", file=sys.stderr)

    # ---- aggregate per discovered leader ------------------------------------
    leader_sites = collections.defaultdict(collections.Counter)
    leader_reads = collections.defaultdict(set)
    leader_tails = collections.defaultdict(set)
    leader_loci = collections.defaultdict(set)
    for s in sites:
        key = s["leader_key"]
        leader_sites[key][(s["chrom"], s["pos"], s["strand"])] += 1
        leader_reads[key].add(s["read"])
        leader_tails[key].add(s["tail"])
        for locus in leaders[key]:
            leader_loci[key].add((locus["chrom"], locus["start"], locus["end"],
                                  locus["strand"]))

    # ---- consolidate fragmented centroids -----------------------------------
    # Tail clustering (SLIDR_CLUSTERTAILS) only merges exact 3'-suffixes.
    # Reads with variable 5' Illumina adapter remnants (e.g. AGATCGGAAGAGCACACGT...)
    # or variable 5' degradation split a single biological leader (e.g. SL1) into
    # hundreds of distinct centroids. Consolidating centroids that share an identical
    # 3'-terminal 15-mer prior to abundance ranking prevents pseudo-leader fragmentation
    # and restores accurate rank ordering.
    candidate_consensus = {}
    candidate_cutoffs = {}
    for key in leader_sites:
        tails_sorted = sorted(leader_tails[key], key=len)
        cons, cut = cluster_consensus(tails_sorted)
        candidate_consensus[key] = cons
        candidate_cutoffs[key] = cut

    if args.consolidation_len > 0:
        kmer_groups = collections.defaultdict(list)
        short_centroids = []

        for key, cons in candidate_consensus.items():
            if len(cons) >= args.consolidation_len:
                kmer = cons[-args.consolidation_len:].upper()
                kmer_groups[kmer].append(key)
            else:
                short_centroids.append((key, cons.upper()))

        # Suffix matching for short centroids (< consolidation_len nt)
        for key, short_seq in short_centroids:
            matched_kmer = None
            best_reads = -1
            for kmer, members in kmer_groups.items():
                if kmer.endswith(short_seq):
                    kmer_reads = sum(len(leader_reads[m]) for m in members)
                    if kmer_reads > best_reads:
                        best_reads = kmer_reads
                        matched_kmer = kmer
            if matched_kmer is not None:
                kmer_groups[matched_kmer].append(key)
            else:
                kmer_groups[short_seq].append(key)

        centroid_groups = list(kmer_groups.values())
    else:
        centroid_groups = [[k] for k in leader_sites]

    reportable = []
    for members in centroid_groups:
        group_reads = set.union(*(leader_reads[m] for m in members))
        group_tails = set.union(*(leader_tails[m] for m in members))
        group_loci = sorted(set.union(*(leader_loci[m] for m in members)))

        group_sites = collections.Counter()
        for m in members:
            for s, c in leader_sites[m].items():
                group_sites[s] += c

        n_reads = len(group_reads)
        strong = {s: c for s, c in group_sites.items() if c >= args.min_support}

        # Mirror the tested rule: a leader must be supported by more than one read
        # and used at more than one independent site, otherwise it is more likely
        # a fragment or a pseudogene than a real leader.
        if n_reads < args.min_support or len(strong) < 2:
            continue

        # Primary anchor: prefer non-adapter centroids, then highest read count, then longest
        members_by_reads = sorted(
            members,
            key=lambda m: (
                not is_adapter_contaminated(candidate_consensus[m]),
                len(leader_reads[m]),
                len(candidate_consensus[m])
            ),
            reverse=True
        )
        best_m = members_by_reads[0]
        top_consensus = candidate_consensus[best_m]
        top_cutoff = candidate_cutoffs[best_m]

        # Extend if a longer non-adapter member with substantial support (>= 10% reads of top)
        # cleanly contains top_consensus as a 3' suffix:
        chosen_consensus = top_consensus
        chosen_cutoff = top_cutoff
        for m in members_by_reads:
            if len(leader_reads[m]) < 0.10 * len(leader_reads[best_m]):
                break
            cons_m = candidate_consensus[m]
            if is_adapter_contaminated(cons_m):
                continue
            if len(cons_m) > len(chosen_consensus) and cons_m.upper().endswith(top_consensus.upper()):
                chosen_consensus = cons_m
                chosen_cutoff = candidate_cutoffs[m]

        family = classify_sl_family(chosen_consensus)

        reportable.append({
            "key": best_m,
            "members": members,
            "consensus": chosen_consensus,
            "cutoff": chosen_cutoff,
            "reads": n_reads,
            "sites": len(strong),
            "site_counts": strong,
            "loci": group_loci,
            "family": family,
        })

    if not reportable:
        print("[slidr_consensus] FATAL: no leader passed the support thresholds "
              f"(>= {args.min_support} reads and >= 2 sites).", file=sys.stderr)
        sys.exit(1)

    # ---- name leaders de novo, most-used first ------------------------------
    reportable.sort(key=lambda r: (-r["sites"], -r["reads"]))
    for idx, rec in enumerate(reportable, 1):
        rec["name"] = f"SL{idx}"

    if rfam_loci:
        annotate_with_rfam(reportable, rfam_loci)

    # ---- outputs -----------------------------------------------------------
    bed_path = os.path.join(args.outdir, f"{args.prefix}_sl_acceptor_sites.bed")
    gff_path = os.path.join(args.outdir, f"{args.prefix}_sl_acceptor_sites.gff3")
    sum_path = os.path.join(args.outdir, f"{args.prefix}_sl_summary.tsv")
    fa_path = os.path.join(args.outdir, f"{args.prefix}_sl_sequences.fa")

    total_sites = 0
    with open(bed_path, "w") as bed, open(gff_path, "w") as gff:
        gff.write("##gff-version 3\n")
        gff.write("# trans-splice acceptor sites and spliced-leader RNA genes "
                  "called de novo by SLIDR\n")
        for rec in reportable:
            name = rec["name"]
            fam = rec["family"]
            rfam_note = f";Rfam_Validation={','.join(rec['rfam_hits'])}" if rec.get("rfam_hits") else ""
            for (chrom, pos, strand), count in sorted(
                    rec["site_counts"].items(), key=lambda kv: (kv[0][0], kv[0][1])):
                total_sites += 1
                # BED6: chromStart is 0-based and chromEnd exclusive, so the
                # junction base at 1-based `pos` becomes [pos-1, pos). The
                # acceptor dinucleotide itself occupies pos-2..pos-1.
                bed.write(f"{chrom}\t{pos - 1}\t{pos}\t{name}:{count}:{fam}\t{count}\t"
                          f"{strand}\n")
                gff.write(f"{chrom}\tSLIDR\ttrans_splice_acceptor_site\t{pos}\t{pos}\t"
                          f"{count}\t{strand}\t.\tID={name}.trans-splice-site-"
                          f"{total_sites};Name={name};Family={fam};Reads={count}{rfam_note}\n")
            for lidx, (chrom, start, end, strand) in enumerate(rec["loci"], 1):
                gff.write(f"{chrom}\tSLIDR\tspliced_leader_RNA\t{start}\t{end}\t.\t"
                          f"{strand}\t.\tID={name}.leader-RNA-{lidx};Name={name};Family={fam}{rfam_note}\n")

    with open(sum_path, "w") as out:
        out.write("sl_name\tconsensus\tlength\treads\tacceptor_sites\t"
                  "leader_RNA_genes\tfamily\trfam_validation\n")
        for rec in reportable:
            rfam_val = ",".join(rec.get("rfam_hits", [])) or "None"
            out.write(f"{rec['name']}\t{rec['consensus']}\t{len(rec['consensus'])}\t"
                      f"{rec['reads']}\t{rec['sites']}\t{len(rec['loci'])}\t{rec['family']}\t{rfam_val}\n")

    with open(fa_path, "w") as out:
        for rec in reportable:
            out.write(f">{rec['name']} reads={rec['reads']} sites={rec['sites']} family={rec['family']}\n"
                      f"{rec['consensus']}\n")

    print(f"[slidr_consensus] Consolidated into {len(reportable)} leaders, {total_sites} acceptor sites "
          f"with >= {args.min_support} supporting reads", file=sys.stderr)
    for rec in reportable[:10]:
        rfam_desc = f" [{','.join(rec.get('rfam_hits', []))}]" if rec.get("rfam_hits") else ""
        print(f"    {rec['name']:<6} ({rec['family']:<18}) reads={rec['reads']:<9} sites={rec['sites']:<7} "
              f"len={len(rec['consensus']):<4} {rec['consensus']}{rfam_desc}", file=sys.stderr)
def load_rfam_sl_loci(path: Optional[str]) -> List[dict]:
    """
    Load curated/predicted SL loci from Rfam cmscan GFF3 (e.g. RF00198 for SL1, RF00199 for SL2).
    """
    if not path or not os.path.exists(path):
        return []
    loci = []
    sl_pat = re.compile(r"RF00198|RF00199|SL1|SL2", re.IGNORECASE)
    with open_maybe_gzip(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            if sl_pat.search(parts[8]) or sl_pat.search(parts[2]):
                try:
                    name_m = re.search(r"Name=([^;]+)", parts[8])
                    acc_m  = re.search(r"Accession=([^;]+)", parts[8])
                    model_name = name_m.group(1) if name_m else parts[2]
                    acc = acc_m.group(1) if acc_m else ""
                    loci.append({
                        "chrom": parts[0],
                        "start": int(parts[3]),
                        "end": int(parts[4]),
                        "strand": parts[6],
                        "model": model_name,
                        "acc": acc
                    })
                except ValueError:
                    continue
    return loci


def annotate_with_rfam(reportable: List[dict], rfam_loci: List[dict], max_overlap_dist: int = 200):
    """
    Cross-reference discovered leaders against Rfam SL loci.
    """
    for rec in reportable:
        rec["rfam_hits"] = []
        for l_chrom, l_start, l_end, l_strand in rec["loci"]:
            for r in rfam_loci:
                if l_chrom == r["chrom"] and l_strand == r["strand"]:
                    if not (l_end < r["start"] - max_overlap_dist or l_start > r["end"] + max_overlap_dist):
                        rec["rfam_hits"].append(f"{r['model']}({r['acc']})")
        rec["rfam_hits"] = sorted(set(rec["rfam_hits"]))


if __name__ == "__main__":
    main()



