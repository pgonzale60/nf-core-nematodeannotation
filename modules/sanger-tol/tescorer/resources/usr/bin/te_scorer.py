#!/usr/bin/env python3
"""
te_scorer.py - Multi-Evidence Tiered Scoring Engine for Eukaryotic Gene Models

Scores predicted gene models on a continuous evidence scale (-100 to +100)
and classifies them into four actionable tiers:
  Tier 1 (Score >= +30) : Confirmed Host Gene
  Tier 2 (Score 0 to 29): Putative Host Gene / Lineage-Specific Orphan
  Tier 3 (Score -40 to -1): Ambiguous / Borderline / TE-derived Host
  Tier 4 (Score < -40)  : Definitive Transposable Element
"""

import os
import argparse
import gzip
from collections import defaultdict

# Specific Pfam / InterPro IDs indicative of Transposable Elements
TE_PFAM_FAMILIES = {
    'PF00078',  # Reverse transcriptase (RNA-dependent DNA polymerase)
    'PF07727',  # Reverse transcriptase (reverse transcriptase domain)
    'PF00665',  # Integrase core domain
    'PF05380',  # Pao retrotransposon peptidase
    'PF01609',  # DDE superfamily endonuclease (Transposase)
    'PF03108',  # Transposase, Tc1-like
    'PF13358',  # DDE superfamily endonuclease / transposase
    'PF02992',  # Transposase family
    'PF03017',  # Transposase, mutator-like
    'PF03732',  # Retrotransposon gag protein
    'PF03733',  # Retrovirus zinc finger-like domain (gag)
    'PF00075',  # RNase H (frequent in retrotransposons)
    'PF01548',  # Maran transposase
    'PF01359',  # Transposase IS4 family
    'PF03400',  # Transposase IS30 family
    'PF13963',  # Transposase
    'PF13843',  # Transposase
    'PF00664',  # Integrase, catalytic core
}

TE_TEXT_KEYWORDS = {
    'reverse transcriptase',
    'retrotransposon',
    'transposase',
    'integrase',
    'gag-pol',
    'pao retrotransposon',
    'mariner transposase',
    'tc1 transposase',
    'retrovirus',
    'endonuclease, reverse transcriptase',
    'rnase h, retroviral'
}

# Negative repeat classes that are genuine TEs in Earl Grey
DEFINED_TE_CLASSES = {
    'line', 'ltr', 'dna', 'sine', 'rc', 'helitron', 'penelope', 'pararetrovirus',
    'ltr/pao', 'ltr/gypsy', 'ltr/copia', 'ltr/bel', 'line/cr1', 'line/l1',
    'line/rte', 'dna/tcmar', 'dna/tcmar-tc1', 'dna/maverick', 'dna/cmc', 'dna/hobo'
}

def parse_args():
    parser = argparse.ArgumentParser(description="Multi-Evidence Tiered TE Scorer")
    parser.add_argument("--gtf", required=True, help="Input consensus GTF (e.g. TSEBRA)")
    parser.add_argument("--fasta", required=True, help="Input predicted protein FASTA")
    parser.add_argument("--iprscan", required=False, help="InterProScan TSV annotations")
    parser.add_argument("--earlgrey", required=False, help="Earl Grey filteredRepeats GFF3")
    parser.add_argument("--satellites", required=False, help="Centroid satellites BED")
    parser.add_argument("--operons", required=False, help="SLOPPR polycistronic operons GFF3")
    parser.add_argument("--orthogroups", required=False, help="OrthoFinder Orthogroups.tsv")
    parser.add_argument("--diamond", required=False, help="Diamond BLASTP TSV against Nematoda reference proteomes")
    parser.add_argument("--out-tsv", default="gene_te_scores.tsv", help="Output scoring table")
    parser.add_argument("--out-host-gtf", default=None, help="Optional output Host GTF (Tiers 1 & 2)")
    parser.add_argument("--out-te-gtf", default=None, help="Optional output TE GTF (Tiers 3 & 4)")
    parser.add_argument("--summary", default="te_scoring_summary.md", help="Summary report")
    return parser.parse_args()

def open_file(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")

def normalize_id(identifier):
    """Normalize transcript / gene identifiers between GTF (anno1.g1.t1) and TSV (g1.t1)."""
    return identifier.replace("anno1.", "").replace("anno2.", "").strip()

def main():
    args = parse_args()

    print("=== Loading Gene Models from GTF ===")
    # Map transcript_id -> {chrom, strand, exons: [(start, end)], gene_id}
    transcripts = {}
    with open_file(args.gtf) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            feature = parts[2]
            if feature not in ("CDS", "exon"):
                continue
            chrom = parts[0]
            start = int(parts[3])
            end = int(parts[4])
            strand = parts[6]
            attr = parts[8]

            tx_id = None
            gene_id = None
            for token in attr.split(";"):
                token = token.strip()
                if token.startswith("transcript_id"):
                    tx_id = token.split('"')[1] if '"' in token else token.split()[1]
                elif token.startswith("gene_id"):
                    gene_id = token.split('"')[1] if '"' in token else token.split()[1]

            if not tx_id:
                continue

            if tx_id not in transcripts:
                transcripts[tx_id] = {
                    'chrom': chrom,
                    'strand': strand,
                    'gene_id': gene_id or tx_id,
                    'cds': [],
                    'exons': []
                }
            if feature == "CDS":
                transcripts[tx_id]['cds'].append((start, end))
            elif feature == "exon":
                transcripts[tx_id]['exons'].append((start, end))

    print(f"Loaded {len(transcripts)} transcripts from GTF.")

    # 1. InterProScan TE and Host Signatures
    tx_te_domains = defaultdict(list)
    tx_host_domains = defaultdict(list)

    if args.iprscan and os.path.exists(args.iprscan):
        print(f"Loading InterProScan annotations from {args.iprscan}...")
        with open_file(args.iprscan) as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) < 6:
                    continue
                tx_raw = parts[0]
                norm_tx = normalize_id(tx_raw)
                db = parts[3]
                sig_acc = parts[4]
                sig_desc = parts[5].lower()
                ipr_acc = parts[11] if len(parts) > 11 else ""
                ipr_desc = parts[12].lower() if len(parts) > 12 else ""

                full_text = f"{sig_acc} {sig_desc} {ipr_acc} {ipr_desc}"

                is_te = False
                if sig_acc in TE_PFAM_FAMILIES or ipr_acc in TE_PFAM_FAMILIES:
                    is_te = True
                elif any(kw in full_text for kw in TE_TEXT_KEYWORDS):
                    is_te = True

                if is_te:
                    tx_te_domains[norm_tx].append(f"{sig_acc}:{parts[5]}")
                else:
                    if ipr_acc != "-" and ipr_acc != "":
                        tx_host_domains[norm_tx].append(f"{ipr_acc}:{ipr_desc}")
        print(f"Found TE domains in {len(tx_te_domains)} transcripts.")

    # 2. OrthoFinder Orthogroups
    # If provided, check if transcript clusters with C. elegans, P. pacificus, or B. malayi
    tx_ortho_hit = {}
    if args.orthogroups and os.path.exists(args.orthogroups):
        print(f"Loading OrthoFinder orthogroups from {args.orthogroups}...")
        with open_file(args.orthogroups) as f:
            header = f.readline().strip().split("\t")
            # header: Orthogroup, sp1, sp2, sp3...
            for line in f:
                parts = line.strip().split("\t")
                og = parts[0]
                has_reference = False
                # Check columns for reference species
                for i in range(1, len(parts)):
                    col_name = header[i].lower()
                    if any(sp in col_name for sp in ['elegans', 'pacificus', 'malayi', 'ref']):
                        if parts[i].strip():
                            has_reference = True
                
                # Assign to all genes in this OG
                for i in range(1, len(parts)):
                    col_name = header[i].lower()
                    if 'arhod' in col_name or 'nxaua' in col_name or 'target' in col_name:
                        for g in parts[i].split(","):
                            g = normalize_id(g.strip())
                            if g:
                                tx_ortho_hit[g] = (og, has_reference)
        print(f"Mapped {len(tx_ortho_hit)} transcripts to Orthogroups.")

    # 3. Diamond BLASTP Fallback / Complementary Homology
    tx_diamond_hit = {}
    if args.diamond and os.path.exists(args.diamond):
        print(f"Loading Diamond BLASTP results from {args.diamond}...")
        with open_file(args.diamond) as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) < 12:
                    continue
                qseqid = normalize_id(parts[0])
                sseqid = parts[1]
                pident = float(parts[2])
                evalue = float(parts[10])
                bitscore = float(parts[11])
                if qseqid not in tx_diamond_hit or bitscore > tx_diamond_hit[qseqid]['bitscore']:
                    tx_diamond_hit[qseqid] = {
                        'target': sseqid,
                        'pident': pident,
                        'evalue': evalue,
                        'bitscore': bitscore
                    }
        print(f"Mapped {len(tx_diamond_hit)} transcripts to Diamond hits.")

    # 4. SLOPPR Operon Demarcation
    tx_operon_status = {}
    if args.operons and os.path.exists(args.operons):
        print(f"Loading SLOPPR operons from {args.operons}...")
        # Check overlaps or explicit gene IDs in GFF3
        operon_intervals = defaultdict(list)
        with open_file(args.operons) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) < 9:
                    continue
                if parts[2] in ("operon", "operonic_gene"):
                    chrom = parts[0]
                    start = int(parts[3])
                    end = int(parts[4])
                    strand = parts[6]
                    operon_intervals[chrom].append((start, end, strand, parts[8]))
        
        # Spatial overlap check against operon coordinates
        for tx_id, info in transcripts.items():
            chrom = info['chrom']
            strand = info['strand']
            if not info['cds']:
                continue
            tx_start = min(s for s, e in info['cds'])
            tx_end = max(e for s, e in info['cds'])
            for op_s, op_e, op_strand, op_attr in operon_intervals.get(chrom, []):
                if strand == op_strand and not (tx_end < op_s or tx_start > op_e):
                    tx_operon_status[tx_id] = "operonic"
                    break

    # 5. Earl Grey Repeat Overlap Calculation
    print("Calculating Earl Grey repeat overlaps...")
    earlgrey_intervals = defaultdict(list)
    if args.earlgrey and os.path.exists(args.earlgrey):
        with open_file(args.earlgrey) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) < 9:
                    continue
                chrom = parts[0]
                rep_type = parts[2].lower()
                start = int(parts[3])
                end = int(parts[4])
                earlgrey_intervals[chrom].append((start, end, rep_type))

    # 6. Satellite Arrays Overlap Calculation
    satellite_intervals = defaultdict(list)
    if args.satellites and os.path.exists(args.satellites):
        print(f"Loading Satellite arrays from {args.satellites}...")
        with open_file(args.satellites) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) < 3:
                    continue
                chrom = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                satellite_intervals[chrom].append((start, end))

    # Evaluate Each Gene Model
    scored_models = []
    tier_counts = {1: 0, 2: 0, 3: 0, 4: 0}

    for tx_id, info in transcripts.items():
        norm_tx = normalize_id(tx_id)
        chrom = info['chrom']
        strand = info['strand']
        gene_id = info['gene_id']
        cds_list = info['cds'] or info['exons']

        total_cds_bp = sum(e - s + 1 for s, e in cds_list)
        tx_start = min(s for s, e in cds_list)
        tx_end = max(e for s, e in cds_list)

        # A. Repeat overlap
        te_overlap_bp = 0
        unknown_overlap_bp = 0
        overlapping_te_classes = set()

        for c_start, c_end in cds_list:
            for r_start, r_end, r_type in earlgrey_intervals.get(chrom, []):
                # Check overlap
                ov_s = max(c_start, r_start)
                ov_e = min(c_end, r_end)
                if ov_s <= ov_e:
                    ov_len = ov_e - ov_s + 1
                    is_defined_te = any(cls in r_type for cls in DEFINED_TE_CLASSES)
                    if is_defined_te:
                        te_overlap_bp += ov_len
                        overlapping_te_classes.add(r_type)
                    else:
                        unknown_overlap_bp += ov_len

        te_pct = (te_overlap_bp / total_cds_bp * 100.0) if total_cds_bp > 0 else 0.0
        unk_pct = (unknown_overlap_bp / total_cds_bp * 100.0) if total_cds_bp > 0 else 0.0

        # B. Satellite overlap
        sat_overlap_bp = 0
        for c_start, c_end in cds_list:
            for s_start, s_end in satellite_intervals.get(chrom, []):
                ov_s = max(c_start, s_start)
                ov_e = min(c_end, s_end)
                if ov_s <= ov_e:
                    sat_overlap_bp += (ov_e - ov_s + 1)
        sat_pct = (sat_overlap_bp / total_cds_bp * 100.0) if total_cds_bp > 0 else 0.0

        # C. Evidence Scoring Computation
        score = 0.0
        reasons = []

        # 1. Homology Score (OrthoFinder / Diamond)
        has_nematode_orthology = False
        if norm_tx in tx_ortho_hit:
            og, has_ref = tx_ortho_hit[norm_tx]
            if has_ref:
                score += 35.0
                has_nematode_orthology = True
                reasons.append(f"OrthoFinder_RefOrtholog({og})[+35]")
            else:
                score += 10.0
                reasons.append(f"OrthoFinder_AuanemaOrtholog({og})[+10]")
        elif norm_tx in tx_diamond_hit:
            hit = tx_diamond_hit[norm_tx]
            if hit['evalue'] < 1e-10:
                score += 30.0
                has_nematode_orthology = True
                reasons.append(f"Diamond_NematodaMatch(E={hit['evalue']:.1e})[+30]")
            elif hit['evalue'] < 1e-5:
                score += 15.0
                reasons.append(f"Diamond_WeakMatch(E={hit['evalue']:.1e})[+15]")

        # 2. Host Functional Domains
        n_host = len(tx_host_domains.get(norm_tx, []))
        if n_host > 0:
            score += min(20.0, 10.0 + n_host * 2.0)
            reasons.append(f"HostDomains({n_host})[+{min(20.0, 10.0 + n_host * 2.0):.0f}]")

        # 3. Operon Context
        if tx_id in tx_operon_status or norm_tx in tx_operon_status:
            score += 15.0
            reasons.append("OperonicGene[+15]")

        # 4. Gene Structure (Exons)
        n_exons = len(cds_list)
        if n_exons >= 3:
            score += 10.0
            reasons.append(f"MultiExon({n_exons})[+10]")
        elif n_exons == 1:
            score -= 10.0
            reasons.append("SingleExon[-10]")

        # 5. TE Functional Domains Penalty
        te_doms = tx_te_domains.get(norm_tx, [])
        if te_doms:
            penalty = 40.0 + min(20.0, len(te_doms) * 5.0)
            score -= penalty
            reasons.append(f"TEDomains({len(te_doms)})[-{penalty:.0f}]")

        # 6. Structural TE Repeat Overlap Penalty
        if te_pct >= 50.0:
            score -= 35.0
            reasons.append(f"HighTEOverlap({te_pct:.1f}%)[-35]")
        elif te_pct >= 20.0:
            score -= 20.0
            reasons.append(f"ModerateTEOverlap({te_pct:.1f}%)[-20]")
        elif te_pct >= 5.0:
            score -= 10.0
            reasons.append(f"LowTEOverlap({te_pct:.1f}%)[-10]")

        # 7. Unknown repeat (slight penalty only if high overlap and NO homology)
        if unk_pct >= 70.0 and not has_nematode_orthology and not n_host:
            score -= 10.0
            reasons.append(f"UnknownRepeatOnly({unk_pct:.1f}%)[-10]")

        # 8. Satellite Overlap Penalty
        if sat_pct >= 20.0:
            score -= 20.0
            reasons.append(f"SatelliteOverlap({sat_pct:.1f}%)[-20]")

        # Determine Tier
        if score >= 30.0:
            tier = 1
            tier_name = "Confirmed_Host"
        elif score >= 0.0:
            tier = 2
            tier_name = "Putative_Host"
        elif score >= -40.0:
            tier = 3
            tier_name = "Ambiguous_TE_Derived"
        else:
            tier = 4
            tier_name = "Definitive_TE"

        tier_counts[tier] += 1
        scored_models.append({
            'transcript_id': tx_id,
            'gene_id': gene_id,
            'chrom': chrom,
            'start': tx_start,
            'end': tx_end,
            'strand': strand,
            'exons': n_exons,
            'cds_bp': total_cds_bp,
            'score': round(score, 1),
            'tier': tier,
            'tier_name': tier_name,
            'te_overlap_pct': round(te_pct, 1),
            'unknown_overlap_pct': round(unk_pct, 1),
            'sat_overlap_pct': round(sat_pct, 1),
            'te_classes': ";".join(sorted(overlapping_te_classes)) if overlapping_te_classes else "None",
            'te_domains': ";".join(te_doms) if te_doms else "None",
            'host_domains_count': n_host,
            'evidence_reasons': " | ".join(reasons)
        })

    # Sort scored models by genomic coordinate
    scored_models.sort(key=lambda m: (m['chrom'], m['start']))

    print("=== Scoring Complete ===")
    print(f"Tier 1 (Confirmed Host)       : {tier_counts[1]}")
    print(f"Tier 2 (Putative Host)        : {tier_counts[2]}")
    print(f"Tier 3 (Ambiguous / Flagged)  : {tier_counts[3]}")
    print(f"Tier 4 (Definitive TE)        : {tier_counts[4]}")

    # Write output scoring TSV
    with open(args.out_tsv, "w") as out:
        headers = [
            "transcript_id", "gene_id", "chrom", "start", "end", "strand",
            "exons", "cds_bp", "score", "tier", "tier_name",
            "te_overlap_pct", "unknown_overlap_pct", "sat_overlap_pct",
            "te_classes", "te_domains", "host_domains_count", "evidence_reasons"
        ]
        out.write("\t".join(headers) + "\n")
        for m in scored_models:
            row = [str(m[h]) for h in headers]
            out.write("\t".join(row) + "\n")

    # Optional Partition GTF into Host vs TE
    if args.out_host_gtf or args.out_te_gtf:
        host_transcripts = {m['transcript_id'] for m in scored_models if m['tier'] in (1, 2)}
        te_transcripts   = {m['transcript_id'] for m in scored_models if m['tier'] in (3, 4)}

        with open_file(args.gtf) as fin:
            f_host = open(args.out_host_gtf, "w") if args.out_host_gtf else None
            f_te = open(args.out_te_gtf, "w") if args.out_te_gtf else None
            for line in fin:
                if line.startswith("#"):
                    if f_host: f_host.write(line)
                    if f_te: f_te.write(line)
                    continue
                parts = line.strip().split("\t")
                if len(parts) < 9:
                    continue
                attr = parts[8]
                tx_id = None
                for token in attr.split(";"):
                    token = token.strip()
                    if token.startswith("transcript_id"):
                        tx_id = token.split('"')[1] if '"' in token else token.split()[1]
                        break
                if tx_id in host_transcripts and f_host:
                    f_host.write(line)
                elif tx_id in te_transcripts and f_te:
                    f_te.write(line)
            if f_host: f_host.close()
            if f_te: f_te.close()

    # Write Markdown Summary
    with open(args.summary, "w") as f_sum:
        total = len(scored_models)
        f_sum.write("# Transposable Element (TE) Tiered Scoring Summary\n\n")
        f_sum.write(f"Evaluated consensus models: **{total:,}**\n\n")
        f_sum.write("| Tier | Classification | Count | Percentage | Action |\n")
        f_sum.write("| :---: | :--- | :---: | :---: | :--- |\n")
        f_sum.write(f"| **Tier 1** | Confirmed Host Gene | {tier_counts[1]:,} | {tier_counts[1]/total*100:.1f}% | Retained in Primary Proteome |\n")
        f_sum.write(f"| **Tier 2** | Putative Host / Orphan | {tier_counts[2]:,} | {tier_counts[2]/total*100:.1f}% | Retained in Primary Proteome |\n")
        f_sum.write(f"| **Tier 3** | Ambiguous / TE-Derived | {tier_counts[3]:,} | {tier_counts[3]/total*100:.1f}% | Flagged / Manual Curation Review |\n")
        f_sum.write(f"| **Tier 4** | Definitive TE Model | {tier_counts[4]:,} | {tier_counts[4]/total*100:.1f}% | Filtered from Host Proteome |\n\n")
        f_sum.write("### Deliverables:\n")
        f_sum.write(f"- Detailed Scoring Sheet: `{args.out_tsv}`\n")
        if args.out_host_gtf:
            f_sum.write(f"- Curated Host GTF: `{args.out_host_gtf}` ({tier_counts[1] + tier_counts[2]:,} models)\n")
        if args.out_te_gtf:
            f_sum.write(f"- Flagged TE Catalog GTF: `{args.out_te_gtf}` ({tier_counts[3] + tier_counts[4]:,} models)\n")

    print(f"Generated {args.out_tsv} and {args.summary}.")

if __name__ == "__main__":
    main()
