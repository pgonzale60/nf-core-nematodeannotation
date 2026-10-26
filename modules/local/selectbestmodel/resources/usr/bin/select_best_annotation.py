#!/usr/bin/env python3
"""
select_best_annotation.py — Run BUSCO on multiple candidate protein sets, compare completeness,
and export the winning annotation and decision table.
"""

import os
import sys
import argparse
import subprocess
import re
import shutil

def run_busco_on_fasta(fasta_path, candidate_name, lineage, lineage_dir, cpus):
    """Run offline protein BUSCO and return complete percentage."""
    if not fasta_path or not os.path.exists(fasta_path):
        return -1.0, "Missing FASTA"

    # Verify fasta has records
    with open(fasta_path, 'r') as f:
        has_records = any(line.startswith('>') for line in f)
    if not has_records:
        return -1.0, "Empty FASTA"

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
        # BUSCO --download_path expects either <download_path>/lineages/<lineage>
        # If user passed a directory containing <lineage> directly instead of lineages/<lineage>,
        # accommodate both by checking and symlinking or pointing as appropriate.
        resolved_download_path = os.path.abspath(lineage_dir)
        direct_lineage_path = os.path.join(resolved_download_path, lineage)
        nested_lineage_path = os.path.join(resolved_download_path, "lineages", lineage)

        if not os.path.exists(nested_lineage_path) and os.path.exists(direct_lineage_path):
            # The download directory directly contains <lineage>, not lineages/<lineage>.
            # Create a local wrapper directory ./busco_lineages_wrapper/lineages/<lineage> via symlink
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
        return -1.0, f"BUSCO execution failed: {e.returncode}"

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
        with open(summary_txt, 'r') as f:
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

    # Clean intermediate files to save inodes
    if os.path.exists(out_prefix):
        for item in os.listdir(out_prefix):
            item_path = os.path.join(out_prefix, item)
            if os.path.isdir(item_path):
                shutil.rmtree(item_path, ignore_errors=True)

    return completeness, summary_str

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

        comp, summary = run_busco_on_fasta(aa, name, args.lineage, args.lineage_dir, args.cpus)
        results.append({
            "name": name,
            "gtf": os.path.abspath(gtf),
            "aa": os.path.abspath(aa),
            "completeness": comp,
            "summary": summary
        })

    if not results:
        print("ERROR: No candidates could be evaluated.")
        sys.exit(1)

    # Sort descending by completeness
    results.sort(key=lambda x: x["completeness"], reverse=True)
    winner = results[0]

    print("\n================== BUSCO TOURNAMENT SUMMARY ==================")
    for r in results:
        star = " ★ (WINNER)" if r["name"] == winner["name"] else ""
        print(f"  {r['name']:<18}: Complete = {r['completeness']:>5.2f}% | {r['summary']}{star}")
    print("==============================================================\n")

    out_gtf = os.path.join(args.outdir, f"{args.prefix}.gtf")
    out_aa = os.path.join(args.outdir, f"{args.prefix}.aa")
    out_tsv = os.path.join(args.outdir, f"{args.prefix}_decision.tsv")

    if os.path.exists(out_gtf) or os.path.islink(out_gtf):
        os.remove(out_gtf)
    if os.path.exists(out_aa) or os.path.islink(out_aa):
        os.remove(out_aa)

    os.symlink(winner['gtf'], out_gtf)
    os.symlink(winner['aa'], out_aa)

    with open(out_tsv, "w") as f:
        f.write("Rank\tCandidate\tBUSCO_Complete_Pct\tSummary\tStatus\tGTF_Path\tAA_Path\n")
        for i, r in enumerate(results, 1):
            status = "SELECTED" if r["name"] == winner["name"] else "REJECTED"
            f.write(f"{i}\t{r['name']}\t{r['completeness']:.2f}\t{r['summary']}\t{status}\t{r['gtf']}\t{r['aa']}\n")

    print(f"Exported winning model to {out_gtf} and {out_aa}")

if __name__ == "__main__":
    main()
