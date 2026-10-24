#!/usr/bin/env python3
"""
sloppr_operons.py — Identify polycistronic operon clusters and operon-demarcated cistrons from trans-splice sites.

Modernized component of SLOPPR (Wenzel et al., BMC Bioinformatics 2021).
Uses trans-splice acceptor sites (specifically SL2 / operon-resolving leaders) and
intercistronic distance constraints to predict operon boundaries and flag internal trans-splicing.
"""

import argparse
import collections
import gzip
import os
import sys
from typing import List

def parse_gff_genes(gff_path: str) -> List[dict]:
    """Parse gene models from GFF3 / GTF."""
    genes = []
    opener = gzip.open if gff_path.endswith(".gz") else open
    with opener(gff_path, "rt") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            feat_type = parts[2].lower()
            if feat_type in ("gene", "mrna", "transcript"):
                chrom = parts[0]
                start = int(parts[3])
                end = int(parts[4])
                strand = parts[6]
                attrs = parts[8]
                gene_id = None
                for attr in attrs.split(";"):
                    if "=" in attr:
                        k, v = attr.strip().split("=", 1)
                        if k.lower() in ("id", "gene_id", "name"):
                            gene_id = v
                            break
                    elif " " in attr.strip():
                        # GTF style
                        k, v = attr.strip().split(" ", 1)
                        if k in ("gene_id", "transcript_id"):
                            gene_id = v.replace('"', '')
                            break
                if not gene_id:
                    gene_id = f"gene_{chrom}_{start}_{end}"
                genes.append({
                    "id": gene_id,
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "strand": strand,
                })
    # Sort genes by chrom, start
    genes.sort(key=lambda g: (g["chrom"], g["start"]))
    return genes


def parse_sl_bed(bed_path: str) -> List[dict]:
    """Parse BED6 trans-splice sites."""
    sites = []
    with open(bed_path, "r") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 6:
                continue
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            name = parts[3]
            score = int(parts[4])
            strand = parts[5]
            sites.append({
                "chrom": chrom,
                "pos": end,  # 1-based acceptor pos
                "name": name,
                "score": score,
                "strand": strand,
                "is_sl2": ("sl2" in name.lower())
            })
    return sites


def main():
    parser = argparse.ArgumentParser(
        description="Predict eukaryotic operons using trans-spliced leader (SL2) acceptor sites."
    )
    parser.add_argument("--bed", required=True, help="Trans-splice acceptor sites BED file (from slidr_consensus.py)")
    parser.add_argument("--gff", default=None, help="Genome annotations GFF3/GTF file (optional)")
    parser.add_argument("--max-intercistronic-dist", type=int, default=1000,
                        help="Maximum intercistronic distance for operon cistrons (default: 1000 bp)")
    parser.add_argument("--prefix", default="aps4", help="Output prefix")
    parser.add_argument("--outdir", default=".", help="Output directory")

    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    sl_sites = parse_sl_bed(args.bed)
    print(f"[sloppr_operons] Loaded {len(sl_sites):,} trans-splice acceptor sites from {args.bed}", file=sys.stderr)

    if not args.gff or not os.path.exists(args.gff):
        print("[sloppr_operons] No gene annotation provided. Grouping trans-splice sites by genomic proximity...", file=sys.stderr)
        # Without gene annotations, cluster proximate SL sites on the same strand
        out_gff = os.path.join(args.outdir, f"{args.prefix}_predicted_operon_clusters.gff3")
        with open(out_gff, "w") as f:
            f.write("##gff-version 3\n")
            # Group by chrom, strand
            by_chrom_strand = collections.defaultdict(list)
            for s in sl_sites:
                by_chrom_strand[(s["chrom"], s["strand"])].append(s)

            cluster_id = 0
            for (chrom, strand), sites in by_chrom_strand.items():
                sites.sort(key=lambda x: x["pos"])
                curr_cluster = []
                for s in sites:
                    if not curr_cluster:
                        curr_cluster.append(s)
                    else:
                        prev = curr_cluster[-1]
                        dist = s["pos"] - prev["pos"]
                        if dist <= args.max_intercistronic_dist:
                            curr_cluster.append(s)
                        else:
                            if len(curr_cluster) > 1 and any(x["is_sl2"] for x in curr_cluster):
                                cluster_id += 1
                                f.write(
                                    f"{chrom}\tSLOPPR\toperon_cluster\t{curr_cluster[0]['pos']}\t{curr_cluster[-1]['pos']}\t"
                                    f".\t{strand}\t.\tID=OP_{cluster_id};Sites={len(curr_cluster)}\n"
                                )
                            curr_cluster = [s]
                if len(curr_cluster) > 1 and any(x["is_sl2"] for x in curr_cluster):
                    cluster_id += 1
                    f.write(
                        f"{chrom}\tSLOPPR\toperon_cluster\t{curr_cluster[0]['pos']}\t{curr_cluster[-1]['pos']}\t"
                        f".\t{strand}\t.\tID=OP_{cluster_id};Sites={len(curr_cluster)}\n"
                    )
        print(f"[sloppr_operons] Emitted {out_gff}", file=sys.stderr)
        return

    # With gene annotations
    genes = parse_gff_genes(args.gff)
    print(f"[sloppr_operons] Loaded {len(genes):,} genes from {args.gff}", file=sys.stderr)

    # Associate SL sites with genes
    # Gene is downstream operonic if its 5' boundary (or 5' UTR / start) receives SL2
    # In C. elegans & nematodes, SL2 marks downstream cistrons in operons.
    operons_gff = os.path.join(args.outdir, f"{args.prefix}_predicted_operons.gff3")
    operons_tsv = os.path.join(args.outdir, f"{args.prefix}_operons_summary.tsv")

    # Index SL sites by (chrom, strand)
    sl_by_chr_strand = collections.defaultdict(list)
    for s in sl_sites:
        sl_by_chr_strand[(s["chrom"], s["strand"])].append(s)

    gene_sl_status = {}
    for g in genes:
        chrom = g["chrom"]
        strand = g["strand"]
        start = g["start"]
        end = g["end"]
        # 5' window: 500 bp upstream to 200 bp downstream of gene start
        if strand == "+":
            w_start = max(1, start - 500)
            w_end = start + 200
        else:
            w_start = end - 200
            w_end = end + 500

        hits = [
            s for s in sl_by_chr_strand.get((chrom, strand), [])
            if w_start <= s["pos"] <= w_end
        ]
        has_sl2 = any(h["is_sl2"] for h in hits)
        has_sl1 = any(not h["is_sl2"] for h in hits)
        gene_sl_status[g["id"]] = {
            "has_sl2": has_sl2,
            "has_sl1": has_sl1,
            "sl_count": len(hits)
        }

    # Identify operon chains
    operon_chains = []
    current_chain = []

    for i in range(len(genes)):
        g = genes[i]
        if not current_chain:
            current_chain.append(g)
            continue

        prev = current_chain[-1]
        same_chrom = (g["chrom"] == prev["chrom"])
        same_strand = (g["strand"] == prev["strand"])
        
        # Intercistronic distance
        if prev["strand"] == "+":
            dist = g["start"] - prev["end"]
        else:
            dist = g["start"] - prev["end"]  # genes sorted by start

        dist = abs(dist)

        # Downstream gene must receive SL2 and be closely linked (<= max distance)
        g_status = gene_sl_status[g["id"]]
        is_downstream_operon = g_status["has_sl2"] and (dist <= args.max_intercistronic_dist)

        if same_chrom and same_strand and is_downstream_operon:
            current_chain.append(g)
        else:
            if len(current_chain) > 1:
                operon_chains.append(current_chain)
            current_chain = [g]

    if len(current_chain) > 1:
        operon_chains.append(current_chain)

    with open(operons_gff, "w") as ogff, open(operons_tsv, "w") as otsv:
        ogff.write("##gff-version 3\n")
        otsv.write("operon_id\tchrom\tstrand\tstart\tend\tnum_cistrons\tcistron_ids\n")

        for idx, chain in enumerate(operon_chains, 1):
            op_id = f"{args.prefix}_OP_{idx}"
            chrom = chain[0]["chrom"]
            strand = chain[0]["strand"]
            op_start = min(g["start"] for g in chain)
            op_end = max(g["end"] for g in chain)
            cistron_ids = ",".join(g["id"] for g in chain)

            ogff.write(
                f"{chrom}\tSLOPPR\toperon\t{op_start}\t{op_end}\t.\t{strand}\t.\t"
                f"ID={op_id};Cistrons={len(chain)};Genes={cistron_ids}\n"
            )
            for c_idx, cistron in enumerate(chain, 1):
                role = "upstream_cistron" if c_idx == 1 else "downstream_cistron"
                ogff.write(
                    f"{chrom}\tSLOPPR\toperonic_gene\t{cistron['start']}\t{cistron['end']}\t.\t{strand}\t.\t"
                    f"ID={cistron['id']};Parent={op_id};Role={role}\n"
                )

            otsv.write(f"{op_id}\t{chrom}\t{strand}\t{op_start}\t{op_end}\t{len(chain)}\t{cistron_ids}\n")

    print(f"[sloppr_operons] Complete! Predicted {len(operon_chains):,} operons across genome.", file=sys.stderr)


if __name__ == "__main__":
    main()
