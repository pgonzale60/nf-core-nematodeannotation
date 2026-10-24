#!/usr/bin/env python3
"""
slidr_tss_exporter.py — Export empirical trans-splicing TSS hints (GFF3).

Consumes the locus-level pool written by SLIDR_CONSENSUS (*_sl_site_pool.tsv), in
which every leader-attributed read is pooled at its genomic acceptor site
irrespective of which tail-cluster centroid it fell into. Thresholding happens
here, on the pooled counts:

    reads        >= --min-reads      independent supporting reads
    n_libraries  >= --min-libraries  independent libraries (1 = no requirement)

No leader identity is assumed. Leader names produced upstream are rank-based
(SL1 = most-used), not biological, so they are deliberately not used to classify
sites. Each hint carries the most frequent 5' tail observed at the site instead.

Outputs:
    {prefix}.tss.gff3         one TSS feature per site, 1-based, strand-aware
    {prefix}.tss_metrics.tsv  counts before/after thresholds
"""

import argparse
import os
import sys


def parse_args():
    p = argparse.ArgumentParser(description="Export empirical trans-splice TSS hints.")
    p.add_argument(
        "--pool", required=True, help="*_sl_site_pool.tsv from SLIDR_CONSENSUS"
    )
    p.add_argument("--prefix", required=True, help="Output file prefix")
    p.add_argument(
        "--min-reads",
        type=int,
        default=3,
        help="Minimum pooled reads per site (default: 3)",
    )
    p.add_argument(
        "--min-libraries",
        type=int,
        default=1,
        help="Minimum independent libraries per site (default: 1)",
    )
    p.add_argument("--outdir", default=".", help="Output directory")
    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    out_gff = os.path.join(args.outdir, f"{args.prefix}.tss.gff3")
    out_tsv = os.path.join(args.outdir, f"{args.prefix}.tss_metrics.tsv")

    counts = {
        "pooled_sites": 0,
        "below_min_reads": 0,
        "below_min_libraries": 0,
        "exported": 0,
        "exported_reads": 0,
    }

    with open(args.pool) as fin, open(out_gff, "w") as fgff:
        fgff.write("##gff-version 3\n")
        fgff.write(
            f"# Empirical trans-splicing TSS hints derived de novo by SLIDR "
            f"for {args.prefix}\n"
        )
        fgff.write(
            f"# thresholds: reads>={args.min_reads} libraries>={args.min_libraries}\n"
        )
        header = fin.readline().rstrip("\n").split("\t")
        ix = {name: i for i, name in enumerate(header)}
        need = (
            "chrom",
            "pos",
            "strand",
            "reads",
            "n_libraries",
            "n_centroids",
            "top_tail",
            "top_tail_reads",
        )
        missing = [c for c in need if c not in ix]
        if missing:
            sys.exit(f"[slidr_tss_exporter] {args.pool} is missing columns {missing}")

        n = 0
        for line in fin:
            p = line.rstrip("\n").split("\t")
            if len(p) < len(header):
                continue
            counts["pooled_sites"] += 1
            reads = int(p[ix["reads"]])
            libs = int(p[ix["n_libraries"]])
            if reads < args.min_reads:
                counts["below_min_reads"] += 1
                continue
            if libs < args.min_libraries:
                counts["below_min_libraries"] += 1
                continue
            n += 1
            counts["exported"] += 1
            counts["exported_reads"] += reads
            pos = p[ix["pos"]]
            attr = (
                f"ID=TSS_{n:06d};reads={reads};libraries={libs};"
                f"centroids={p[ix['n_centroids']]};"
                f"top_tail={p[ix['top_tail']]};"
                f"top_tail_reads={p[ix['top_tail_reads']]};src=trans_splicing"
            )
            fgff.write(
                f"{p[ix['chrom']]}\tSLIDR\tTSS\t{pos}\t{pos}\t{reads}\t"
                f"{p[ix['strand']]}\t.\t{attr}\n"
            )

    with open(out_tsv, "w") as ftsv:
        ftsv.write("metric\tcount\n")
        ftsv.writelines(f"{k}\t{v}\n" for k, v in counts.items())

    print(
        f"[slidr_tss_exporter] {counts['exported']:,} TSS hints exported from "
        f"{counts['pooled_sites']:,} pooled sites "
        f"(reads>={args.min_reads}, libraries>={args.min_libraries}) -> {out_gff}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
