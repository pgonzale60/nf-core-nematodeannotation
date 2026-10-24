#!/usr/bin/env python3
"""
filter_trf_satellites.py — Filter TRF output for high-confidence satellite repeats.

Criteria:
- Unit length (period size) >= min_period (default: 50 bp)
- Copy number >= min_copies (default: 100.0)
- Sequence identity (percent match) >= min_identity (default: 85%)

Outputs:
- <prefix>_filtered_units.fa: Consensus monomer sequences for clustering by VSEARCH
- <prefix>_trf_summary.tsv: Tabular details of qualifying repeat loci
"""

import argparse
import os
import sys

def parse_trf_dat(dat_file):
    """Parse TRF .dat output format."""
    current_seq = "unknown"
    repeats = []
    
    with open(dat_file, "r") as f:
        for line in f:
            line_str = line.strip()
            if not line_str:
                continue
            if line_str.startswith("Sequence:"):
                parts_seq = line_str.split()
                if len(parts_seq) > 1:
                    current_seq = parts_seq[1]
                continue
            if line_str.startswith("Parameters:"):
                continue
            
            parts = line_str.split()
            # TRF tabular repeat line has at least 14 columns
            # 0: start, 1: end, 2: period, 3: copies, 4: consensus_size, 5: pct_matches, 6: pct_indels,
            # 7: score, 8: A%, 9: C%, 10: G%, 11: T%, 12: entropy, 13: consensus_seq
            if len(parts) >= 14:
                try:
                    start = int(parts[0])
                    end = int(parts[1])
                    period = int(parts[2])
                    copies = float(parts[3])
                    consensus_size = int(parts[4])
                    pct_matches = int(parts[5])
                    pct_indels = int(parts[6])
                    score = int(parts[7])
                    consensus = parts[13].upper()
                    
                    repeats.append({
                        "seq_id": current_seq,
                        "start": start,
                        "end": end,
                        "period": period,
                        "copies": copies,
                        "consensus_size": consensus_size,
                        "pct_matches": pct_matches,
                        "pct_indels": pct_indels,
                        "score": score,
                        "consensus": consensus
                    })
                except ValueError:
                    continue
    return repeats

def main():
    parser = argparse.ArgumentParser(description="Filter TRF outputs for satellite repeat consensus harvesting.")
    parser.add_argument("--dat", nargs="+", required=True, help="One or more TRF .dat output files")
    parser.add_argument("--min-period", type=int, default=50, help="Minimum period size in bp (default: 50)")
    parser.add_argument("--min-copies", type=float, default=100.0, help="Minimum copy number (default: 100.0)")
    parser.add_argument("--min-identity", type=int, default=85, help="Minimum percent match identity (default: 85)")
    parser.add_argument("--prefix", default="satellites", help="Output prefix")
    parser.add_argument("--outdir", default=".", help="Output directory")

    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    out_fa = os.path.join(args.outdir, f"{args.prefix}_filtered_units.fa")
    out_tsv = os.path.join(args.outdir, f"{args.prefix}_trf_summary.tsv")

    total_parsed = 0
    passed_repeats = []

    for dat_path in args.dat:
        if os.path.isfile(dat_path):
            reps = parse_trf_dat(dat_path)
            total_parsed += len(reps)
            for r in reps:
                if (r["period"] >= args.min_period and 
                    r["copies"] >= args.min_copies and 
                    r["pct_matches"] >= args.min_identity):
                    passed_repeats.append(r)

    print(f"[filter_trf_satellites] Evaluated {total_parsed:,} total TRF repeats across {len(args.dat)} file(s).", file=sys.stderr)
    print(f"[filter_trf_satellites] {len(passed_repeats):,} passed criteria (period>={args.min_period}bp, copies>={args.min_copies}, id>={args.min_identity}%).", file=sys.stderr)

    with open(out_fa, "w") as fa_out, open(out_tsv, "w") as tsv_out:
        tsv_out.write("seq_id\tstart\tend\tperiod\tcopies\tpct_matches\tscore\tconsensus_len\tconsensus_seq\n")
        for idx, r in enumerate(passed_repeats, 1):
            unit_id = f"{args.prefix}_sat_{idx}_{r['seq_id']}_{r['start']}_{r['end']}_p{r['period']}_c{int(r['copies'])}"
            fa_out.write(f">{unit_id}\n{r['consensus']}\n")
            tsv_out.write(f"{r['seq_id']}\t{r['start']}\t{r['end']}\t{r['period']}\t{r['copies']}\t{r['pct_matches']}\t{r['score']}\t{len(r['consensus'])}\t{r['consensus']}\n")

    print(f"[filter_trf_satellites] Wrote {len(passed_repeats):,} consensus units to {out_fa}", file=sys.stderr)
    print(f"[filter_trf_satellites] Summary table saved to {out_tsv}", file=sys.stderr)

if __name__ == "__main__":
    main()
