#!/usr/bin/env python3
"""
satellite_profiler.py - Tandem Repeat & Satellite Architecture Profiler

Purpose:
  Transforms raw BLASTN alignments of clustered tandem repeat centroids (from TRF+VSEARCH)
  into a fully annotated, family-classified satellite catalog:
  1. Canonical Circular Phasing: Normalizes monomer rotations to the lexicographical
     minimum so circularly permuted monomers share identical phase.
  2. Family Abundance Ranking: Measures total genome coverage per centroid and assigns
     ranked, clean family identifiers (sat_fam01, sat_fam02, ...).
  3. Array Resolution: Clusters adjacent/overlapping hits of the same family into discrete
     Satellite Arrays (calculating array span, copy number, and purity %).
  4. Dual GFF3 & Table Emission: Outputs both array-level and monomer-level GFF3 annotations,
     a curated monomer consensus FASTA library, and a comprehensive summary table.
"""

import argparse
import gzip
import os
import sys
from collections import defaultdict


def open_file(path, mode="rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def canonical_rotation(seq: str) -> str:
    """
    Returns the lexicographically smallest circular shift of seq.
    Runs Booth's algorithm-equivalent via suffix comparison in O(n) time.
    """
    s = seq.upper()
    n = len(s)
    ss = s + s
    min_rot = s
    min_idx = 0
    # Booth's Algorithm
    f = [-1] * (2 * n)
    k = 0
    for j in range(1, 2 * n):
        sj = ss[j]
        i = f[j - k - 1]
        while i != -1 and sj != ss[k + i + 1]:
            if sj < ss[k + i + 1]:
                k = j - i - 1
            i = f[i]
        if i == -1 and sj != ss[k]:
            if sj < ss[k]:
                k = j
            f[j - k] = -1
        else:
            f[j - k] = i + 1
    return ss[k : k + n]


def parse_centroids_fasta(fasta_path: str):
    """
    Parses centroids FASTA and returns dict: {centroid_id: sequence}
    """
    centroids = {}
    current_id = None
    seq_chunks = []
    with open_file(fasta_path, "rt") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id:
                    centroids[current_id] = "".join(seq_chunks)
                current_id = line[1:].split()[0]
                seq_chunks = []
            else:
                seq_chunks.append(line)
        if current_id:
            centroids[current_id] = "".join(seq_chunks)
    return centroids


def parse_blast_bed(bed_path: str):
    """
    Parses BLAST-to-BED intervals:
    Format: chrom, start, end, centroid_id, pident, strand
    """
    records = []
    with open_file(bed_path, "rt") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 6:
                continue
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            centroid_id = parts[3]
            try:
                pident = float(parts[4])
            except ValueError:
                pident = 100.0
            strand = parts[5] if parts[5] in ("+", "-") else "+"
            records.append({
                "chrom": chrom,
                "start": start,
                "end": end,
                "centroid_id": centroid_id,
                "pident": pident,
                "strand": strand,
                "length": end - start
            })
    return records


def cluster_into_arrays(monomers, max_gap=100):
    """
    Merges sorted monomer hits of the same family into contiguous Satellite Arrays.
    """
    if not monomers:
        return []

    arrays = []
    current_array = [monomers[0]]

    for m in monomers[1:]:
        prev = current_array[-1]
        if m["chrom"] == prev["chrom"] and (m["start"] - prev["end"]) <= max_gap:
            current_array.append(m)
        else:
            arrays.append(current_array)
            current_array = [m]
    if current_array:
        arrays.append(current_array)

    resolved_arrays = []
    for arr in arrays:
        chrom = arr[0]["chrom"]
        arr_start = arr[0]["start"]
        arr_end = arr[-1]["end"]
        arr_len = arr_end - arr_start
        copy_count = len(arr)
        total_pident = sum(m["pident"] for m in arr)
        mean_pident = round(total_pident / copy_count, 1)

        # Monomer widths
        widths = [m["end"] - m["start"] for m in arr]
        widths.sort()
        median_width = widths[len(widths) // 2]

        resolved_arrays.append({
            "chrom": chrom,
            "start": arr_start,
            "end": arr_end,
            "length": arr_len,
            "copy_count": copy_count,
            "mean_pident": mean_pident,
            "median_monomer_len": median_width,
            "monomers": arr
        })
    return resolved_arrays


def main():
    parser = argparse.ArgumentParser(
        description="Profile, classify, and annotate satellite DNA arrays from TRF/VSEARCH/BLASTN."
    )
    parser.add_argument("--bed", required=True, help="Sorted BLASTN-derived BED file (chrom, start, end, centroid_id, pident, strand)")
    parser.add_argument("--centroids", required=True, help="FASTA of clustered centroids from VSEARCH")
    parser.add_argument("--prefix", default="assembly", help="Output prefix")
    parser.add_argument("--max-gap", type=int, default=100, help="Max gap (bp) to merge monomers into an array (default: 100)")
    parser.add_argument("--outdir", default=".", help="Output directory")

    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # 1. Load centroids and monomers
    centroids = parse_centroids_fasta(args.centroids)
    monomers = parse_blast_bed(args.bed)

    if not monomers:
        sys.stderr.write("No satellite monomers found in input BED.\n")
        # Touch empty deliverables
        for suffix in ("_satellite_arrays.gff3", "_satellite_monomers.gff3", "_satellite_families.tsv", "_satellite_monomers.fasta"):
            with open(os.path.join(args.outdir, f"{args.prefix}{suffix}"), "w") as f:
                pass
        return 0

    # 2. Measure coverage per centroid and rank families
    centroid_coverage = defaultdict(int)
    centroid_counts = defaultdict(int)
    for m in monomers:
        centroid_coverage[m["centroid_id"]] += m["length"]
        centroid_counts[m["centroid_id"]] += 1

    ranked_centroids = sorted(
        centroid_coverage.keys(), key=lambda c: centroid_coverage[c], reverse=True
    )

    family_map = {}
    family_metadata = {}
    for idx, c_id in enumerate(ranked_centroids, 1):
        fam_name = f"sat_fam{idx:02d}"
        family_map[c_id] = fam_name

        raw_seq = centroids.get(c_id, "")
        canon_seq = canonical_rotation(raw_seq) if raw_seq else ""
        period = len(canon_seq)
        gc_pct = round(((canon_seq.count("G") + canon_seq.count("C")) / period * 100), 1) if period > 0 else 0.0

        family_metadata[fam_name] = {
            "centroid_id": c_id,
            "canonical_seq": canon_seq,
            "period": period,
            "gc_pct": gc_pct,
            "total_bp": centroid_coverage[c_id],
            "copy_count": centroid_counts[c_id],
        }

    # Assign family name to monomers
    for m in monomers:
        m["family"] = family_map.get(m["centroid_id"], "sat_unclassified")

    # 3. Group monomers by family and cluster into discrete arrays
    monomers_by_fam = defaultdict(list)
    for m in monomers:
        monomers_by_fam[m["family"]].append(m)

    all_arrays = []
    for fam, f_monomers in monomers_by_fam.items():
        # Ensure chromosome order then coordinate order
        f_monomers.sort(key=lambda x: (x["chrom"], x["start"]))
        arrays = cluster_into_arrays(f_monomers, max_gap=args.max_gap)
        for a in arrays:
            a["family"] = fam
            all_arrays.append(a)

    # Sort arrays globally by genomic position
    all_arrays.sort(key=lambda x: (x["chrom"], x["start"]))

    # 4. Write Deliverable: satellite_families.tsv
    summary_tsv_path = os.path.join(args.outdir, f"{args.prefix}_satellite_families.tsv")
    with open(summary_tsv_path, "w") as f:
        f.write("family_id\tcentroid_id\tmonomer_length_bp\tgc_percent\ttotal_coverage_bp\tcopy_count\tarray_count\n")
        fam_array_counts = defaultdict(int)
        for a in all_arrays:
            fam_array_counts[a["family"]] += 1

        for fam_name in sorted(family_metadata.keys()):
            meta = family_metadata[fam_name]
            arr_count = fam_array_counts[fam_name]
            f.write(
                f"{fam_name}\t{meta['centroid_id']}\t{meta['period']}\t{meta['gc_pct']}\t"
                f"{meta['total_bp']}\t{meta['copy_count']}\t{arr_count}\n"
            )

    # 5. Write Deliverable: satellite_monomers.fasta (Phase-normalized consensus library)
    library_fasta_path = os.path.join(args.outdir, f"{args.prefix}_satellite_monomers.fasta")
    with open(library_fasta_path, "w") as f:
        for fam_name in sorted(family_metadata.keys()):
            meta = family_metadata[fam_name]
            f.write(
                f">{fam_name} period={meta['period']}bp gc={meta['gc_pct']}% total_cov={meta['total_bp']}bp copies={meta['copy_count']}\n"
            )
            seq = meta["canonical_seq"]
            # 60 chars per line
            for i in range(0, len(seq), 60):
                f.write(seq[i : i + 60] + "\n")

    # 6. Write Deliverable: satellite_arrays.gff3
    arrays_gff3_path = os.path.join(args.outdir, f"{args.prefix}_satellite_arrays.gff3")
    with open(arrays_gff3_path, "w") as f:
        f.write("##gff-version 3\n")
        for i, a in enumerate(all_arrays, 1):
            attr = (
                f"ID=array_{i:04d};Family={a['family']};"
                f"Copies={a['copy_count']};MedianUnitLen={a['median_monomer_len']};"
                f"MeanIdentity={a['mean_pident']}"
            )
            # 1-based start for GFF3
            gff_start = a["start"] + 1
            gff_end = a["end"]
            f.write(
                f"{a['chrom']}\tTRF_Satellite\tsatellite_array\t{gff_start}\t{gff_end}\t.\t.\t.\t{attr}\n"
            )

    # 7. Write Deliverable: satellite_monomers.gff3
    monomers_gff3_path = os.path.join(args.outdir, f"{args.prefix}_satellite_monomers.gff3")
    with open(monomers_gff3_path, "w") as f:
        f.write("##gff-version 3\n")
        monomer_counter = 1
        for arr_idx, a in enumerate(all_arrays, 1):
            for m in a["monomers"]:
                attr = (
                    f"ID=monomer_{monomer_counter:06d};Parent=array_{arr_idx:04d};"
                    f"Family={m['family']};Identity={m['pident']}"
                )
                gff_start = m["start"] + 1
                gff_end = m["end"]
                f.write(
                    f"{m['chrom']}\tTRF_Satellite\tsatellite_monomer\t{gff_start}\t{gff_end}\t{m['pident']}\t{m['strand']}\t.\t{attr}\n"
                )
                monomer_counter += 1

    print(f"Satellite profiling complete: {len(family_metadata)} families, {len(all_arrays)} arrays, {len(monomers)} monomers.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
