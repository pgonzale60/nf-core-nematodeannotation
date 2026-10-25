#!/usr/bin/env python3
"""
filter_structural_tandem_arrays.py: Extract high-density tandem arrays of structural RNAs
(Pol I rRNA, 5S rRNA, spliceosomal snRNAs) into a BED mask to shield TE discovery
engines (e.g., Earl Grey / RepeatModeler) from runaway compute, while leaving dispersed,
solitary tRNAs, 7SL, and single retrotransposed snRNA loci unmasked to preserve genuine biology.
"""

import sys
import argparse
import re

# Structural classes and their clustering thresholds:
# - Major rDNA (SSU 18S, LSU 28S, 5.8S): unit size ~7-12kb including IGS/ITS. Bridge up to 10kb.
# - 5S rRNA: unit size ~120bp + ~1kb spacer. Bridge up to 2kb.
# - Spliceosomal snRNAs (U1, U2, U4, U5, U6, U3): unit size ~100-200bp + short spacers. Bridge up to 1kb.

def parse_gff(gff_path):
    rdna = []
    five_s = []
    snrna = []

    with open(gff_path, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9 or parts[2] != 'ncRNA':
                continue

            chrom = parts[0]
            start = int(parts[3]) - 1 # 0-based
            end = int(parts[4])

            # Parse attributes
            attrs = {}
            for item in parts[8].split(';'):
                if '=' in item:
                    k, v = item.split('=', 1)
                    attrs[k] = v

            name = attrs.get('Name', '')

            if any(k in name for k in ['SSU_rRNA', 'LSU_rRNA', '5_8S_rRNA']):
                rdna.append((chrom, start, end, name))
            elif name == '5S_rRNA':
                five_s.append((chrom, start, end, name))
            elif re.match(r'^U[1-6]$', name) or name == 'U3':
                snrna.append((chrom, start, end, name))

    return rdna, five_s, snrna

def cluster_features(features, max_gap, min_copies, min_span):
    if not features:
        return []

    features.sort(key=lambda x: (x[0], x[1]))
    clusters = []

    for chrom, s, e, name in features:
        if not clusters or clusters[-1][0] != chrom or (s - clusters[-1][2]) > max_gap:
            clusters.append([chrom, s, e, 1, [name]])
        else:
            clusters[-1][2] = max(clusters[-1][2], e)
            clusters[-1][3] += 1
            clusters[-1][4].append(name)

    # Filter by copy count or minimum array span
    return [c for c in clusters if c[3] >= min_copies or (c[2] - c[1]) >= min_span]

def main():
    parser = argparse.ArgumentParser(description="Extract structural RNA tandem arrays for pre-EarlGrey masking")
    parser.add_argument("--gff", required=True, help="Input ncRNA GFF3 from Infernal cmscan")
    parser.add_argument("--out", required=True, help="Output BED file of tandem arrays to hard-mask")
    args = parser.parse_args()

    rdna, five_s, snrna = parse_gff(args.gff)

    # 1. rDNA arrays: Bridge across ITS/IGS up to 10kb; require >= 3 copies or >= 10kb span
    rdna_clusters = cluster_features(rdna, max_gap=10000, min_copies=3, min_span=10000)

    # 2. 5S rRNA: Bridge up to 2kb; require >= 5 tandem copies
    five_s_clusters = cluster_features(five_s, max_gap=2000, min_copies=5, min_span=2000)

    # 3. Spliceosomal snRNAs: Bridge up to 1kb; require >= 5 tandem copies
    snrna_clusters = cluster_features(snrna, max_gap=1000, min_copies=5, min_span=2000)

    all_clusters = rdna_clusters + five_s_clusters + snrna_clusters
    all_clusters.sort(key=lambda x: (x[0], x[1]))

    # Merge any overlapping intervals across structural classes
    merged = []
    for c in all_clusters:
        chrom, s, e = c[0], c[1], c[2]
        if not merged or merged[-1][0] != chrom or s > merged[-1][2]:
            merged.append([chrom, s, e])
        else:
            merged[-1][2] = max(merged[-1][2], e)

    with open(args.out, 'w', encoding='utf-8') as out_f:
        for chrom, s, e in merged:
            out_f.write(f"{chrom}\t{s}\t{e}\n")

if __name__ == '__main__':
    main()
