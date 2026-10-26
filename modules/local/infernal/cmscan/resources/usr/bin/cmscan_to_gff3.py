#!/usr/bin/env python3
"""
cmscan_to_gff3.py — Convert Infernal cmscan --tblout to GFF3 and exclusion BED.
"""

import argparse
import os
import sys

def parse_cmscan_tbl(tbl_file, max_evalue=0.01):
    records = []
    with open(tbl_file, "r") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.strip().split()
            if len(parts) < 17:
                continue
            target_name = parts[0]
            target_acc  = parts[1]
            query_name  = parts[2]
            # Coordinates in cmscan tblout: seq_from (idx 7), seq_to (idx 8), strand (idx 9)
            try:
                seq_from = int(parts[7])
                seq_to   = int(parts[8])
                strand   = parts[9]
                score    = float(parts[14])
                evalue   = float(parts[15])
                inc      = parts[16] # '!' for significant, '?' for inclusion threshold
            except ValueError:
                continue

            if evalue > max_evalue:
                continue

            start = min(seq_from, seq_to)
            end   = max(seq_from, seq_to)
            desc  = " ".join(parts[17:]) if len(parts) > 17 else target_name

            records.append({
                "chrom": query_name,
                "start": start,
                "end": end,
                "strand": strand,
                "score": score,
                "evalue": evalue,
                "target": target_name,
                "acc": target_acc,
                "desc": desc
            })
    return records

def main():
    parser = argparse.ArgumentParser(description="Convert Infernal cmscan tblout to GFF3 and BED.")
    parser.add_argument("--tbl", required=True, help="Infernal cmscan --tblout file")
    parser.add_argument("--prefix", default="aps4_ncrna", help="Output prefix")
    parser.add_argument("--max-evalue", type=float, default=0.01, help="Max E-value threshold")
    parser.add_argument("--outdir", default=".", help="Output directory")

    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    records = parse_cmscan_tbl(args.tbl, max_evalue=args.max_evalue)
    records.sort(key=lambda r: (r["chrom"], r["start"]))

    gff3_path = os.path.join(args.outdir, f"{args.prefix}.gff3")
    bed_path  = os.path.join(args.outdir, f"{args.prefix}_exclusion.bed")

    with open(gff3_path, "w") as gff, open(bed_path, "w") as bed:
        gff.write("##gff-version 3\n")
        for idx, r in enumerate(records, 1):
            feat_id = f"ncRNA_{args.prefix}_{idx}"
            gff.write(
                f"{r['chrom']}\tInfernal\tncRNA\t{r['start']}\t{r['end']}\t{r['score']}\t{r['strand']}\t.\t"
                f"ID={feat_id};Name={r['target']};Accession={r['acc']};Evalue={r['evalue']};Note={r['desc']}\n"
            )
            # 0-based BED for exclusion mask
            bed.write(f"{r['chrom']}\t{r['start'] - 1}\t{r['end']}\t{feat_id}\t{int(r['score'])}\t{r['strand']}\n")

    print(f"[cmscan_to_gff3] Emitted {len(records):,} ncRNA features to {gff3_path} and {bed_path}", file=sys.stderr)

if __name__ == "__main__":
    main()
