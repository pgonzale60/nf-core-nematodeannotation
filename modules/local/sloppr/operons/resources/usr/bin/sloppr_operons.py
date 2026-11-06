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
import re
import sys
from typing import Dict, List, Optional


def classify_sl_family(seq: str) -> str:
    """
    Classify spliced leader sequence into biological families based on 3' conserved core motif:
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


def load_leader_family_map(fasta_path: Optional[str] = None, summary_path: Optional[str] = None) -> Dict[str, str]:
    """Load mapping of leader name -> biological family from FASTA or summary TSV."""
    fam_map = {}
    if summary_path and os.path.exists(summary_path):
        opener = gzip.open if summary_path.endswith(".gz") else open
        with opener(summary_path, "rt") as f:
            header = f.readline().rstrip("\n").split("\t")
            name_idx = 0
            seq_idx = 1 if len(header) > 1 else None
            fam_idx = header.index("family") if "family" in header else None
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) <= name_idx:
                    continue
                name = parts[name_idx]
                if fam_idx is not None and len(parts) > fam_idx:
                    fam_map[name] = parts[fam_idx]
                elif seq_idx is not None and len(parts) > seq_idx:
                    fam_map[name] = classify_sl_family(parts[seq_idx])
    elif fasta_path and os.path.exists(fasta_path):
        opener = gzip.open if fasta_path.endswith(".gz") else open
        with opener(fasta_path, "rt") as f:
            curr_name = None
            curr_seq = []
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if curr_name and curr_seq:
                        seq_str = "".join(curr_seq)
                        fam_map[curr_name] = classify_sl_family(seq_str)
                    header_line = line[1:].split()[0]
                    curr_name = header_line
                    fam_match = re.search(r"family=([A-Za-z0-9_]+)", line)
                    if fam_match:
                        fam_map[curr_name] = fam_match.group(1)
                    curr_seq = []
                else:
                    curr_seq.append(line)
            if curr_name and curr_seq and curr_name not in fam_map:
                fam_map[curr_name] = classify_sl_family("".join(curr_seq))
    return fam_map


def parse_gff_genes(gff_path: str) -> List[dict]:
    """
    Parse deduplicated canonical gene models from GFF3 / GTF.
    Aggregates multi-isoform transcripts by gene_id to prevent single genes from chaining into false self-operons.
    """
    gene_records = collections.defaultdict(lambda: {
        "chrom": None,
        "strand": None,
        "starts": [],
        "ends": [],
    })

    opener = gzip.open if gff_path.endswith(".gz") else open

    # Check if explicit 'gene' features exist in the file
    has_gene_features = False
    with opener(gff_path, "rt") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 3 and parts[2].lower() == "gene":
                has_gene_features = True
                break

    target_features = {"gene"} if has_gene_features else {"mrna", "transcript", "exon", "cds"}

    with opener(gff_path, "rt") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            feat_type = parts[2].lower()
            if feat_type not in target_features:
                continue

            chrom = parts[0]
            start = int(parts[3])
            end = int(parts[4])
            strand = parts[6]
            attrs = parts[8]

            gene_id = None
            if has_gene_features and feat_type == "gene":
                for attr in attrs.split(";"):
                    if "=" in attr:
                        k, v = attr.strip().split("=", 1)
                        if k.lower() in ("id", "gene_id", "name"):
                            gene_id = v
                            break
                    elif " " in attr.strip():
                        k, v = attr.strip().split(" ", 1)
                        if k in ("gene_id", "gene_name"):
                            gene_id = v.replace('"', '')
                            break
            else:
                for attr in attrs.split(";"):
                    if "=" in attr:
                        k, v = attr.strip().split("=", 1)
                        if k.lower() in ("parent", "gene_id", "gene_name"):
                            gene_id = v
                            break
                    elif " " in attr.strip():
                        k, v = attr.strip().split(" ", 1)
                        if k in ("gene_id", "gene_name"):
                            gene_id = v.replace('"', '')
                            break
                if not gene_id:
                    for attr in attrs.split(";"):
                        if "=" in attr:
                            k, v = attr.strip().split("=", 1)
                            if k.lower() in ("id", "transcript_id"):
                                gene_id = v
                                break
                        elif " " in attr.strip():
                            k, v = attr.strip().split(" ", 1)
                            if k == "transcript_id":
                                gene_id = v.replace('"', '')
                                break

            if not gene_id:
                gene_id = f"gene_{chrom}_{start}_{end}"

            rec = gene_records[gene_id]
            rec["chrom"] = chrom
            rec["strand"] = strand
            rec["starts"].append(start)
            rec["ends"].append(end)

    genes = []
    for gid, data in gene_records.items():
        genes.append({
            "id": gid,
            "chrom": data["chrom"],
            "start": min(data["starts"]),
            "end": max(data["ends"]),
            "strand": data["strand"],
        })

    # Sort genes by chrom, start
    genes.sort(key=lambda g: (g["chrom"], g["start"]))
    return genes


def parse_sl_bed(bed_path: str, leader_family_map: Optional[Dict[str, str]] = None) -> List[dict]:
    """Parse BED6 trans-splice sites with biological leader classification."""
    sites = []
    if leader_family_map is None:
        leader_family_map = {}
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
            raw_name = parts[3]
            score = int(parts[4])
            strand = parts[5]

            name_parts = raw_name.split(":")
            lead_id = name_parts[0]

            family = "Unknown"
            if len(name_parts) >= 3 and name_parts[2] in (
                "SL1_canonical", "SL2_canonical_like", "SL_operon_variant",
                "PolyA_PolyT_artifact", "Other_minor"
            ):
                family = name_parts[2]
            elif lead_id in leader_family_map:
                family = leader_family_map[lead_id]

            # In nematodes, downstream operon cistrons are trans-spliced by SL2 or operon-specific leaders
            is_operon_sl = (family in ("SL2_canonical_like", "SL_operon_variant"))

            # Fallback only when no family annotation or sequence table is provided:
            # Match canonical SL2 identifier (e.g. SL2 or SL2_xxx), never matching SL20..SL29!
            if not is_operon_sl and family == "Unknown":
                if re.fullmatch(r"SL2(_.*)?", lead_id, re.IGNORECASE):
                    is_operon_sl = True

            sites.append({
                "chrom": chrom,
                "pos": end,  # 1-based acceptor pos
                "name": lead_id,
                "score": score,
                "strand": strand,
                "family": family,
                "is_sl2": is_operon_sl,
            })
    return sites


def main():
    parser = argparse.ArgumentParser(
        description="Predict eukaryotic operons using trans-spliced leader (SL2) acceptor sites."
    )
    parser.add_argument("--bed", required=True, help="Trans-splice acceptor sites BED file (from slidr_consensus.py)")
    parser.add_argument("--gff", default=None, help="Genome annotations GFF3/GTF file (optional)")
    parser.add_argument("--sl-fasta", default=None, help="Optional SL consensus FASTA file for sequence motif matching")
    parser.add_argument("--sl-summary", default=None, help="Optional SL summary TSV file for family lookup")
    parser.add_argument("--max-intercistronic-dist", type=int, default=1000,
                        help="Maximum intercistronic distance for operon cistrons (default: 1000 bp)")
    parser.add_argument("--prefix", default="aps4", help="Output prefix")
    parser.add_argument("--outdir", default=".", help="Output directory")

    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    fam_map = load_leader_family_map(args.sl_fasta, args.sl_summary)
    # If not provided, inspect directory of bed file for *_sl_sequences.fa or *_sl_summary.tsv
    if not fam_map and args.bed:
        bed_dir = os.path.dirname(os.path.abspath(args.bed))
        for fname in os.listdir(bed_dir):
            if fname.endswith("_sl_sequences.fa") or fname.endswith("_sl_summary.tsv"):
                fpath = os.path.join(bed_dir, fname)
                fam_map.update(load_leader_family_map(
                    fasta_path=fpath if fname.endswith(".fa") else None,
                    summary_path=fpath if fname.endswith(".tsv") else None
                ))
    if fam_map:
        print(f"[sloppr_operons] Loaded {len(fam_map)} leader family classifications.", file=sys.stderr)

    sl_sites = parse_sl_bed(args.bed, fam_map)
    print(f"[sloppr_operons] Loaded {len(sl_sites):,} trans-splice acceptor sites from {args.bed}", file=sys.stderr)

    if not args.gff or not os.path.exists(args.gff):
        print("[sloppr_operons] No gene annotation provided. Grouping trans-splice sites by genomic proximity...", file=sys.stderr)
        out_gff = os.path.join(args.outdir, f"{args.prefix}_predicted_operon_clusters.gff3")
        with open(out_gff, "w") as f:
            f.write("##gff-version 3\n")
            by_chrom_strand = collections.defaultdict(list)
            for s in sl_sites:
                by_chrom_strand[(s["chrom"], s["strand"])].append(s)

            cluster_id = 0
            for (chrom, strand), sites in by_chrom_strand.items():
                if strand == "+":
                    sites.sort(key=lambda x: x["pos"])
                else:
                    sites.sort(key=lambda x: -x["pos"])
                curr_cluster = []
                for s in sites:
                    if not curr_cluster:
                        curr_cluster.append(s)
                    else:
                        prev = curr_cluster[-1]
                        dist = abs(s["pos"] - prev["pos"])
                        if dist <= args.max_intercistronic_dist:
                            curr_cluster.append(s)
                        else:
                            if len(curr_cluster) > 1 and any(x["is_sl2"] for x in curr_cluster):
                                cluster_id += 1
                                c_start = min(x["pos"] for x in curr_cluster)
                                c_end = max(x["pos"] for x in curr_cluster)
                                f.write(
                                    f"{chrom}\tSLOPPR\toperon_cluster\t{c_start}\t{c_end}\t"
                                    f".\t{strand}\t.\tID=OP_{cluster_id};Sites={len(curr_cluster)}\n"
                                )
                            curr_cluster = [s]
                if len(curr_cluster) > 1 and any(x["is_sl2"] for x in curr_cluster):
                    cluster_id += 1
                    c_start = min(x["pos"] for x in curr_cluster)
                    c_end = max(x["pos"] for x in curr_cluster)
                    f.write(
                        f"{chrom}\tSLOPPR\toperon_cluster\t{c_start}\t{c_end}\t"
                        f".\t{strand}\t.\tID=OP_{cluster_id};Sites={len(curr_cluster)}\n"
                    )
        print(f"[sloppr_operons] Emitted {out_gff}", file=sys.stderr)
        return

    # With gene annotations
    genes = parse_gff_genes(args.gff)
    print(f"[sloppr_operons] Loaded {len(genes):,} deduplicated canonical genes from {args.gff}", file=sys.stderr)

    # Associate SL sites with genes
    # Gene is downstream operonic if its 5' boundary receives SL2 / operon-resolving leader
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
        # 5' window: 500 bp upstream to 250 bp downstream of gene start
        if strand == "+":
            w_start = max(1, start - 500)
            w_end = start + 250
        else:
            w_start = end - 250
            w_end = end + 500

        hits = [
            s for s in sl_by_chr_strand.get((chrom, strand), [])
            if w_start <= s["pos"] <= w_end
        ]
        has_sl2 = any(h["is_sl2"] for h in hits)
        has_sl1 = any(h["family"] == "SL1_canonical" for h in hits)
        gene_sl_status[g["id"]] = {
            "has_sl2": has_sl2,
            "has_sl1": has_sl1,
            "sl_count": len(hits)
        }

    # Identify operon chains per chromosome and strand
    operon_chains = []
    chromosomes = sorted(list(set(g["chrom"] for g in genes)))

    for chrom in chromosomes:
        # (+) strand: sorted ascending by start (transcription left to right)
        plus_genes = [g for g in genes if g["chrom"] == chrom and g["strand"] == "+"]
        plus_genes.sort(key=lambda x: x["start"])

        chain = []
        for g in plus_genes:
            if not chain:
                chain.append(g)
                continue

            prev = chain[-1]
            dist = g["start"] - prev["end"]
            g_status = gene_sl_status[g["id"]]
            is_downstream_operon = g_status["has_sl2"] and (0 <= dist <= args.max_intercistronic_dist)

            if is_downstream_operon:
                chain.append(g)
            else:
                if len(chain) > 1:
                    operon_chains.append(chain)
                chain = [g]

        if len(chain) > 1:
            operon_chains.append(chain)

        # (-) strand: sorted descending by end (transcription right to left)
        minus_genes = [g for g in genes if g["chrom"] == chrom and g["strand"] == "-"]
        minus_genes.sort(key=lambda x: -x["end"])

        chain = []
        for g in minus_genes:
            if not chain:
                chain.append(g)
                continue

            prev = chain[-1]
            dist = prev["start"] - g["end"]
            g_status = gene_sl_status[g["id"]]
            is_downstream_operon = g_status["has_sl2"] and (0 <= dist <= args.max_intercistronic_dist)

            if is_downstream_operon:
                chain.append(g)
            else:
                if len(chain) > 1:
                    operon_chains.append(chain)
                chain = [g]

        if len(chain) > 1:
            operon_chains.append(chain)

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
