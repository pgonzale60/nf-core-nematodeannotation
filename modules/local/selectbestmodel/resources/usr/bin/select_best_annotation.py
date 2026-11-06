#!/usr/bin/env python3
"""
select_best_annotation.py — Run BUSCO on multiple candidate protein sets, compare completeness,
filter hits by longest protein isoform per gene to eliminate artificial duplication inflation,
and export the winning representative annotation, BUSCO full table, and decision table.
"""

import os
import sys
import argparse
import subprocess
import re
import shutil
from collections import defaultdict


def parse_fasta_lengths(fasta_path):
    """Parse protein FASTA file and return a dictionary of {sequence_id: length}."""
    lengths = {}
    if not fasta_path or not os.path.exists(fasta_path):
        return lengths
    current_id = None
    current_len = 0
    with open(fasta_path, 'r', encoding='utf-8', errors='replace') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_id is not None:
                    lengths[current_id] = current_len
                current_id = line[1:].split()[0]
                current_len = 0
            else:
                current_len += len(line)
        if current_id is not None:
            lengths[current_id] = current_len
    return lengths


def parse_gtf_gene_transcripts(gtf_path):
    """
    Parse candidate GTF (or GFF3) to build:
    1. transcript_to_gene: mapping from transcript_id / protein_id -> gene_id
    2. gene_to_transcripts: mapping from gene_id -> set of transcript_ids
    3. cds_lengths: mapping from transcript_id -> total CDS length in bp
    """
    transcript_to_gene = {}
    gene_to_transcripts = defaultdict(set)
    cds_lengths = defaultdict(int)

    if not gtf_path or not os.path.exists(gtf_path):
        return transcript_to_gene, gene_to_transcripts, cds_lengths

    gene_id_re = re.compile(r'gene_id\s+["\']?([^"\';]+)["\']?')
    transcript_id_re = re.compile(r'transcript_id\s+["\']?([^"\';]+)["\']?')
    protein_id_re = re.compile(r'protein_id\s+["\']?([^"\';]+)["\']?')

    # GFF3 fallback patterns
    id_re = re.compile(r'(?:^|;)ID=([^;]+)')
    parent_re = re.compile(r'(?:^|;)Parent=([^;]+)')

    with open(gtf_path, 'r', encoding='utf-8', errors='replace') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\r\n').split('\t')
            if len(parts) < 9:
                continue
            feature = parts[2]
            try:
                start = int(parts[3])
                end = int(parts[4])
                length = abs(end - start) + 1
            except ValueError:
                length = 0
            attrs = parts[8]

            gene_id = None
            transcript_id = None
            protein_id = None

            m_gene = gene_id_re.search(attrs)
            if m_gene:
                gene_id = m_gene.group(1).strip()
            m_trans = transcript_id_re.search(attrs)
            if m_trans:
                transcript_id = m_trans.group(1).strip()
            m_prot = protein_id_re.search(attrs)
            if m_prot:
                protein_id = m_prot.group(1).strip()

            # GFF3 fallback
            if not gene_id or not transcript_id:
                m_id = id_re.search(attrs)
                m_parent = parent_re.search(attrs)
                if feature in ('mRNA', 'transcript'):
                    if m_id and not transcript_id:
                        transcript_id = m_id.group(1).strip()
                    if m_parent and not gene_id:
                        gene_id = m_parent.group(1).strip()
                elif feature == 'CDS':
                    if m_parent and not transcript_id:
                        transcript_id = m_parent.group(1).strip()

            if gene_id and transcript_id:
                transcript_to_gene[transcript_id] = gene_id
                gene_to_transcripts[gene_id].add(transcript_id)
                if protein_id:
                    transcript_to_gene[protein_id] = gene_id
                    gene_to_transcripts[gene_id].add(protein_id)

            if feature == 'CDS' and transcript_id and length > 0:
                cds_lengths[transcript_id] += length
                if protein_id:
                    cds_lengths[protein_id] += length

    return transcript_to_gene, gene_to_transcripts, cds_lengths


def get_longest_isoform_per_gene(gtf_path, fasta_path):
    """
    Determine the longest protein isoform for each gene.
    Returns:
    - longest_isoforms: set of transcript/protein IDs designated as the longest isoform
    - transcript_to_gene: mapping of transcript/protein ID -> gene_id
    """
    transcript_to_gene, gene_to_transcripts, cds_lengths = parse_gtf_gene_transcripts(gtf_path)
    fasta_lengths = parse_fasta_lengths(fasta_path)

    longest_isoforms = set()

    for gene_id, transcripts in gene_to_transcripts.items():
        if not transcripts:
            continue
        # Sort transcripts by length descending, with transcript ID tie-breaker
        sorted_transcripts = sorted(
            transcripts,
            key=lambda t: (fasta_lengths.get(t, cds_lengths.get(t, 0)), t),
            reverse=True
        )
        longest_t = sorted_transcripts[0]
        longest_isoforms.add(longest_t)

    # For sequences in FASTA that weren't in GTF (orphans), treat as their own gene and keep as representative
    for seq_id in fasta_lengths:
        if seq_id not in transcript_to_gene:
            transcript_to_gene[seq_id] = seq_id
            longest_isoforms.add(seq_id)

    return longest_isoforms, transcript_to_gene


def calculate_longest_isoform_busco(full_table_path, longest_isoforms, transcript_to_gene):
    """
    Parse BUSCO full_table.tsv, filter hits to retain only the longest protein isoform
    per gene, and recalculate deduplicated BUSCO completeness metrics (C, S, D, F, M).

    Returns:
    (completeness_pct, summary_str, details_dict)
    """
    if not full_table_path or not os.path.exists(full_table_path):
        return -1.0, "Missing full_table.tsv", {}

    hits_by_busco = defaultdict(list)
    all_busco_ids = set()
    total_searched_from_header = None

    with open(full_table_path, 'r', encoding='utf-8', errors='replace') as f:
        for line in f:
            line_str = line.strip()
            if not line_str:
                continue
            if line_str.startswith('#'):
                m_tot = re.search(r'Number of Total BUSCO groups searched:\s*(\d+)', line_str)
                if m_tot:
                    total_searched_from_header = int(m_tot.group(1))
                continue

            parts = line_str.split('\t')
            if len(parts) < 2:
                continue

            busco_id = parts[0].strip()
            status = parts[1].strip()
            seq_id = parts[2].strip() if len(parts) > 2 else ""
            score = float(parts[3].strip()) if len(parts) > 3 and parts[3].strip() not in ("", "N/A", "-") else 0.0
            length = float(parts[4].strip()) if len(parts) > 4 and parts[4].strip() not in ("", "N/A", "-") else 0.0

            all_busco_ids.add(busco_id)
            hits_by_busco[busco_id].append({
                "status": status,
                "seq_id": seq_id,
                "score": score,
                "length": length
            })

    total_buscos = total_searched_from_header if total_searched_from_header else len(all_busco_ids)
    if total_buscos == 0:
        return -1.0, "No BUSCO markers found in full table", {}

    single_count = 0
    duplicated_count = 0
    fragmented_count = 0
    missing_count = 0

    for busco_id in all_busco_ids:
        hits = hits_by_busco.get(busco_id, [])

        complete_genes = set()
        fragmented_genes = set()

        for h in hits:
            status = h["status"]
            seq_id = h["seq_id"]

            if not seq_id:
                continue

            is_longest = (seq_id in longest_isoforms)
            if not is_longest:
                norm_seq = re.sub(r'^(transcript|protein):', '', seq_id)
                if norm_seq in longest_isoforms:
                    is_longest = True
                    seq_id = norm_seq

            if not is_longest:
                continue

            gene_id = transcript_to_gene.get(seq_id, seq_id)

            if status in ("Complete", "Duplicated"):
                complete_genes.add(gene_id)
            elif status == "Fragmented":
                fragmented_genes.add(gene_id)

        if len(complete_genes) == 1:
            single_count += 1
        elif len(complete_genes) > 1:
            duplicated_count += 1
        elif len(fragmented_genes) >= 1:
            fragmented_count += 1
        else:
            missing_count += 1

    complete_count = single_count + duplicated_count
    c_pct = (complete_count / total_buscos) * 100.0
    s_pct = (single_count / total_buscos) * 100.0
    d_pct = (duplicated_count / total_buscos) * 100.0
    f_pct = (fragmented_count / total_buscos) * 100.0
    m_pct = (missing_count / total_buscos) * 100.0

    summary_str = f"C:{c_pct:.1f}%[S:{s_pct:.1f}%,D:{d_pct:.1f}%],F:{f_pct:.1f}%,M:{m_pct:.1f}%"
    details = {
        "C": complete_count,
        "S": single_count,
        "D": duplicated_count,
        "F": fragmented_count,
        "M": missing_count,
        "total": total_buscos,
        "C_pct": c_pct,
        "S_pct": s_pct,
        "D_pct": d_pct,
        "F_pct": f_pct,
        "M_pct": m_pct
    }
    return c_pct, summary_str, details


def run_busco_on_fasta(fasta_path, candidate_name, lineage, lineage_dir, cpus):
    """Run offline protein BUSCO and return (complete_percentage, summary_str, saved_full_table_path)."""
    if not fasta_path or not os.path.exists(fasta_path):
        return -1.0, "Missing FASTA", None

    # Verify fasta has records
    with open(fasta_path, 'r', encoding='utf-8', errors='replace') as f:
        has_records = any(line.startswith('>') for line in f)
    if not has_records:
        return -1.0, "Empty FASTA", None

    out_prefix = f"busco_{candidate_name}"
    cmd = [
        "busco",
        "--cpu", str(cpus),
        "--in", fasta_path,
        "--out", out_prefix,
        "--mode", "proteins",
        "--lineage_dataset", lineage,
        "--offline"
    ]
    if lineage_dir:
        resolved_download_path = os.path.abspath(lineage_dir)
        direct_lineage_path = os.path.join(resolved_download_path, lineage)
        nested_lineage_path = os.path.join(resolved_download_path, "lineages", lineage)

        if not os.path.exists(nested_lineage_path) and os.path.exists(direct_lineage_path):
            wrapper_dir = os.path.abspath("busco_lineages_wrapper")
            wrapper_lineages = os.path.join(wrapper_dir, "lineages")
            os.makedirs(wrapper_lineages, exist_ok=True)
            target_link = os.path.join(wrapper_lineages, lineage)
            if not os.path.exists(target_link):
                try:
                    os.symlink(direct_lineage_path, target_link)
                except OSError:
                    pass
            resolved_download_path = wrapper_dir

        cmd.extend(["--download_path", resolved_download_path])

    print(f"=== Running BUSCO on {candidate_name} ({fasta_path}) ===")
    try:
        subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    except subprocess.CalledProcessError as e:
        print(f"WARNING: BUSCO failed for {candidate_name}: {e.stderr[:300]}")
        return -1.0, f"BUSCO execution failed: {e.returncode}", None

    # Parse short summary
    summary_txt = None
    for root, _, files in os.walk(out_prefix):
        for file in files:
            if file.startswith("short_summary") and file.endswith(".txt"):
                summary_txt = os.path.join(root, file)
                break
        if summary_txt:
            break

    completeness = -1.0
    summary_str = "No summary found"
    if summary_txt and os.path.exists(summary_txt):
        with open(summary_txt, 'r', encoding='utf-8', errors='replace') as f:
            content = f.read()
            m = re.search(r'C:(\d+\.?\d*)%\[S:(\d+\.?\d*)%,D:(\d+\.?\d*)%\],F:(\d+\.?\d*)%,M:(\d+\.?\d*)%', content)
            if m:
                completeness = float(m.group(1))
                summary_str = m.group(0)
            else:
                m_simple = re.search(r'C:(\d+\.?\d*)%', content)
                if m_simple:
                    completeness = float(m_simple.group(1))
                    summary_str = m_simple.group(0)

    # Locate and preserve full_table.tsv before intermediate cleanup
    full_table_tsv = None
    for root, _, files in os.walk(out_prefix):
        for file in files:
            if file == "full_table.tsv" or (file.startswith("full_table") and file.endswith(".tsv")):
                full_table_tsv = os.path.join(root, file)
                break
        if full_table_tsv:
            break

    saved_full_table = f"busco_{candidate_name}_full_table.tsv"
    if full_table_tsv and os.path.exists(full_table_tsv):
        shutil.copyfile(full_table_tsv, saved_full_table)
    else:
        saved_full_table = None

    # Clean intermediate files to save inodes
    if os.path.exists(out_prefix):
        for item in os.listdir(out_prefix):
            item_path = os.path.join(out_prefix, item)
            if os.path.isdir(item_path):
                shutil.rmtree(item_path, ignore_errors=True)

    return completeness, summary_str, saved_full_table


def main():
    parser = argparse.ArgumentParser(description="Evaluate multiple gene predictions via BUSCO and select best.")
    parser.add_argument("--lineage", default="chromadorea_odb12")
    parser.add_argument("--lineage-dir", default="")
    parser.add_argument("--cpus", type=int, default=4)
    parser.add_argument("--prefix", default="representative_annotation")
    parser.add_argument("--outdir", default=".")

    # Candidate pairs: --candidate Name gtf aa
    parser.add_argument("--candidate", action="append", nargs=3, metavar=("NAME", "GTF", "AA"),
                        help="Candidate model: --candidate NAME path/to/gtf path/to/aa")

    args = parser.parse_args()

    if not args.candidate:
        print("ERROR: No candidates provided.")
        sys.exit(1)

    results = []
    for name, gtf, aa in args.candidate:
        if not gtf or gtf == "null" or not os.path.exists(gtf):
            print(f"Skipping {name}: GTF not available ({gtf})")
            continue
        if not aa or aa == "null" or not os.path.exists(aa):
            print(f"Skipping {name}: AA not available ({aa})")
            continue

        comp, summary, cand_full_table = run_busco_on_fasta(aa, name, args.lineage, args.lineage_dir, args.cpus)

        # Compute deduplicated longest-isoform BUSCO score
        li_comp = -1.0
        li_summary = "N/A"
        if cand_full_table and os.path.exists(cand_full_table):
            longest_isoforms, transcript_to_gene = get_longest_isoform_per_gene(gtf, aa)
            li_comp, li_summary, _ = calculate_longest_isoform_busco(cand_full_table, longest_isoforms, transcript_to_gene)

        results.append({
            "name": name,
            "gtf": os.path.abspath(gtf),
            "aa": os.path.abspath(aa),
            "completeness": comp,
            "summary": summary,
            "longest_isoform_completeness": li_comp,
            "longest_isoform_summary": li_summary,
            "full_table": cand_full_table
        })

    if not results:
        print("ERROR: No candidates could be evaluated.")
        sys.exit(1)

    # Sort descending by deduplicated completeness if available, breaking ties with raw completeness
    results.sort(
        key=lambda x: (
            x["longest_isoform_completeness"] if x["longest_isoform_completeness"] >= 0 else x["completeness"],
            x["completeness"]
        ),
        reverse=True
    )
    winner = results[0]

    print("\n================== BUSCO TOURNAMENT SUMMARY ==================")
    for r in results:
        star = " ★ (WINNER)" if r["name"] == winner["name"] else ""
        li_c = r.get("longest_isoform_completeness", -1.0)
        li_sum = r.get("longest_isoform_summary", "N/A")
        print(f"  {r['name']:<18}: Raw C = {r['completeness']:>5.2f}% ({r['summary']}) | Longest-iso C = {li_c:>5.2f}% ({li_sum}){star}")
    print("==============================================================\n")

    out_gtf = os.path.join(args.outdir, f"{args.prefix}.gtf")
    out_aa = os.path.join(args.outdir, f"{args.prefix}.aa")
    out_tsv = os.path.join(args.outdir, f"{args.prefix}_decision.tsv")
    out_busco_table = os.path.join(args.outdir, f"{args.prefix}_busco_full_table.tsv")

    if os.path.exists(out_gtf) or os.path.islink(out_gtf):
        os.remove(out_gtf)
    if os.path.exists(out_aa) or os.path.islink(out_aa):
        os.remove(out_aa)
    if os.path.exists(out_busco_table) or os.path.islink(out_busco_table):
        os.remove(out_busco_table)

    shutil.copy2(winner['gtf'], out_gtf)
    shutil.copy2(winner['aa'], out_aa)

    # Retain BUSCO full_table.tsv for the winning candidate annotation
    if winner.get("full_table") and os.path.exists(winner["full_table"]):
        shutil.copyfile(winner["full_table"], out_busco_table)
    else:
        with open(out_busco_table, "w") as f:
            f.write("# Busco id\tStatus\tSequence\tScore\tLength\n")

    # Record both raw and deduplicated (longest-isoform) BUSCO scores in decision table
    with open(out_tsv, "w") as f:
        f.write("Rank\tCandidate\tBUSCO_Complete_Pct\tSummary\tLongest_Isoform_BUSCO_Complete_Pct\tLongest_Isoform_Summary\tStatus\tGTF_Path\tAA_Path\n")
        for i, r in enumerate(results, 1):
            status = "SELECTED" if r["name"] == winner["name"] else "REJECTED"
            raw_c = r["completeness"]
            raw_sum = r["summary"]
            li_c = r.get("longest_isoform_completeness", -1.0)
            li_sum = r.get("longest_isoform_summary", "N/A")
            f.write(f"{i}\t{r['name']}\t{raw_c:.2f}\t{raw_sum}\t{li_c:.2f}\t{li_sum}\t{status}\t{r['gtf']}\t{r['aa']}\n")

    print(f"Exported winning model to {out_gtf} and {out_aa}")
    print(f"Exported representative BUSCO full table to {out_busco_table}")
    print(f"Exported decision table to {out_tsv}")


if __name__ == "__main__":
    main()
