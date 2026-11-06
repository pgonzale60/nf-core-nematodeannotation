#!/usr/bin/env python3
"""Authoritative data contracts and schema validation for nf-core/nematodeannotation.

Validates:
  1. `assets/schemas/inter_module_tsv_contracts.json` exists, is valid JSON, and
     contains expected top-level schema metadata and contract entries.
  2. Live Python module code matches contract column specifications and status enums:
     - `select_best_annotation.py` -> `decision.tsv` columns and statuses
     - `te_scorer.py` -> `gene_te_scores.tsv` columns and statuses
     - `sloppr_operons.py` -> `operons_summary.tsv` columns
     - `satellite_profiler.py` -> `satellite_families.tsv` columns
     - `filter_trf_satellites.py` -> `trf_summary.tsv` columns
     - `slidr_cluster_tails.py` -> `cluster_summary.tsv` columns
     - `slidr_consensus.py` -> `sl_summary.tsv` and `sl_site_pool.tsv` columns
     - `slidr_extract_tails.py` -> `acceptor_sites.tsv.gz` columns and `tail_stats.json` keys
     - `slidr_filter_sl_rna.py` -> `sl_rna_candidates.tsv` columns
     - `slidr_tss_exporter.py` -> `tss_metrics.tsv` columns
  3. Module `stub:` blocks across `modules/local/**/main.nf` create expected output files
     and emit valid stub headers matching contracts.
  4. Both `scripts/extract_code_contracts.py --check` and `scripts/stamp_module_shas.py --check`
     pass with exit code 0.
"""

from __future__ import annotations

import json
from pathlib import Path
import re
import subprocess
import sys


REPO_ROOT = Path(__file__).resolve().parent.parent
TSV_CONTRACTS_PATH = REPO_ROOT / "assets" / "schemas" / "inter_module_tsv_contracts.json"


def extract_tsv_header_from_file(file_path: Path, pattern: str) -> list[str]:
    """Extract column list from write/write_header pattern in source code."""
    assert file_path.is_file(), f"Source file does not exist: {file_path}"
    content = file_path.read_text(encoding="utf-8")
    match = re.search(pattern, content, re.DOTALL)
    assert match, f"Could not match pattern '{pattern}' in {file_path}"
    raw_header = match.group(1).strip()
    return [col.strip() for col in raw_header.split(r"\t") if col.strip()]


def extract_list_literal(file_path: Path, pattern: str) -> list[str]:
    """Extract list of string literals from python source."""
    assert file_path.is_file(), f"Source file does not exist: {file_path}"
    content = file_path.read_text(encoding="utf-8")
    match = re.search(pattern, content, re.DOTALL)
    assert match, f"Could not match list literal '{pattern}' in {file_path}"
    items = re.findall(r"['\"]([^'\"]+)['\"]", match.group(1))
    assert items, f"No items extracted from {file_path} for pattern {pattern}"
    return items


def verify_schemas_structure() -> dict:
    """Validate existence, JSON syntax, and key requirements for schema files."""
    assert TSV_CONTRACTS_PATH.is_file(), f"Missing schema: {TSV_CONTRACTS_PATH}"

    try:
        tsv_contracts = json.loads(TSV_CONTRACTS_PATH.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise AssertionError(f"Malformed JSON in {TSV_CONTRACTS_PATH}: {exc}") from exc

    assert tsv_contracts.get("schema_version") == "1.0.0", (
        f"Unexpected schema_version in {TSV_CONTRACTS_PATH}: {tsv_contracts.get('schema_version')}"
    )
    assert "description" in tsv_contracts, f"Missing 'description' in {TSV_CONTRACTS_PATH}"
    contracts = tsv_contracts.get("contracts", {})
    assert isinstance(contracts, dict), f"'contracts' must be an object in {TSV_CONTRACTS_PATH}"

    expected_contracts = {
        "tournament_decision.tsv",
        "decision.tsv",
        "busco_full_table.tsv",
        "gene_te_scores.tsv",
        "te_summary.tsv",
        "operons_summary.tsv",
        "operons.tsv",
        "satellite_families.tsv",
        "satellite_summary.tsv",
        "trf_summary.tsv",
        "cluster_summary.tsv",
        "readinfo.tsv.gz",
        "sl_summary.tsv",
        "sl_site_pool.tsv",
        "sl_acceptor_sites.bed",
        "acceptor_sites.tsv.gz",
        "sl_rna_candidates.tsv",
        "tss_metrics.tsv",
        "tail_stats.json",
    }
    missing_contracts = expected_contracts - set(contracts.keys())
    assert not missing_contracts, f"Missing contracts in {TSV_CONTRACTS_PATH}: {missing_contracts}"

    print("✓ Schema file exists, is valid JSON, and meets structural specifications.")
    return tsv_contracts


def verify_live_code_contracts(tsv_contracts: dict) -> None:
    """Verify live Python module code against contract column specifications."""
    contracts = tsv_contracts["contracts"]

    # 1. SELECTBESTMODEL (tournament_decision.tsv and busco_full_table.tsv)
    sbm_path = (
        REPO_ROOT
        / "modules"
        / "local"
        / "selectbestmodel"
        / "resources"
        / "usr"
        / "bin"
        / "select_best_annotation.py"
    )
    sbm_content = sbm_path.read_text(encoding="utf-8")

    # tournament_decision.tsv / decision.tsv
    m_dec = re.search(r'out_tsv[^\n]*\n\s*f\.write\("([^\n]+?)\\n"\)', sbm_content)
    assert m_dec, f"Could not match decision header in {sbm_path}"
    sbm_dec_cols = [c.strip() for c in m_dec.group(1).split(r"\t") if c.strip()]
    expected_dec = contracts["tournament_decision.tsv"]["columns"]
    assert sbm_dec_cols == expected_dec, (
        f"Schema drift in {sbm_path} (decision): {sbm_dec_cols} != {expected_dec}"
    )

    # busco_full_table.tsv
    m_busco = re.search(r'out_busco_table[^\n]*\n\s*f\.write\("([^\n]+?)\\n"\)', sbm_content)
    assert m_busco, f"Could not match busco header in {sbm_path}"
    sbm_busco_cols = [c.strip() for c in m_busco.group(1).split(r"\t") if c.strip()]
    expected_busco = contracts["busco_full_table.tsv"]["columns"]
    assert sbm_busco_cols == expected_busco, (
        f"Schema drift in {sbm_path} (busco): {sbm_busco_cols} != {expected_busco}"
    )

    assert '"SELECTED"' in sbm_content and '"REJECTED"' in sbm_content, (
        f"Expected statuses 'SELECTED' and 'REJECTED' not found in {sbm_path}"
    )
    print("✓ SELECTBESTMODEL (tournament_decision.tsv & busco_full_table.tsv) match live code.")

    # 2. TESCORER (gene_te_scores.tsv)
    tesc_path = (
        REPO_ROOT
        / "modules"
        / "local"
        / "tescorer"
        / "resources"
        / "usr"
        / "bin"
        / "te_scorer.py"
    )
    tesc_cols = extract_list_literal(tesc_path, r"headers\s*=\s*\[(.*?)\]\s*\n\s*out\.write")
    expected_tesc = contracts["gene_te_scores.tsv"]["columns"]
    assert tesc_cols == expected_tesc, (
        f"Schema drift in {tesc_path}: code columns {tesc_cols} != contract columns {expected_tesc}"
    )
    tesc_content = tesc_path.read_text(encoding="utf-8")
    for st in contracts["gene_te_scores.tsv"]["allowed_statuses"]:
        assert f'"{st}"' in tesc_content, f"Expected tier name '{st}' not found in {tesc_path}"
    print("✓ TESCORER (gene_te_scores.tsv) columns and statuses match live code.")

    # 3. SLOPPR_OPERONS (operons_summary.tsv)
    slop_path = (
        REPO_ROOT
        / "modules"
        / "local"
        / "sloppr"
        / "operons"
        / "resources"
        / "usr"
        / "bin"
        / "sloppr_operons.py"
    )
    slop_cols = extract_tsv_header_from_file(slop_path, r'otsv\.write\("([^\n]+?\\t[^\n]+?)\\n"\)')
    expected_slop = contracts["operons_summary.tsv"]["columns"]
    assert slop_cols == expected_slop, (
        f"Schema drift in {slop_path}: code columns {slop_cols} != contract columns {expected_slop}"
    )
    print("✓ SLOPPR_OPERONS (operons_summary.tsv) columns match live code.")

    # 4. SATELLITEPROFILER (satellite_families.tsv)
    satp_path = (
        REPO_ROOT
        / "modules"
        / "local"
        / "satelliteprofiler"
        / "resources"
        / "usr"
        / "bin"
        / "satellite_profiler.py"
    )
    satp_cols = extract_tsv_header_from_file(satp_path, r'f\.write\("([^\n]+?\\t[^\n]+?)\\n"\)')
    expected_satp = contracts["satellite_families.tsv"]["columns"]
    assert satp_cols == expected_satp, (
        f"Schema drift in {satp_path}: code columns {satp_cols} != contract columns {expected_satp}"
    )
    print("✓ SATELLITEPROFILER (satellite_families.tsv) columns match live code.")

    # 5. SATELLITEFILTER (trf_summary.tsv)
    satf_path = (
        REPO_ROOT
        / "modules"
        / "local"
        / "satellitefilter"
        / "resources"
        / "usr"
        / "bin"
        / "filter_trf_satellites.py"
    )
    satf_cols = extract_tsv_header_from_file(satf_path, r'tsv_out\.write\("([^\n]+?\\t[^\n]+?)\\n"\)')
    expected_satf = contracts["trf_summary.tsv"]["columns"]
    assert satf_cols == expected_satf, (
        f"Schema drift in {satf_path}: code columns {satf_cols} != contract columns {expected_satf}"
    )
    print("✓ SATELLITEFILTER (trf_summary.tsv) columns match live code.")

    # 6. SLIDR_CLUSTERTAILS (cluster_summary.tsv)
    slc_path = (
        REPO_ROOT
        / "modules"
        / "local"
        / "slidr"
        / "clustertails"
        / "resources"
        / "usr"
        / "bin"
        / "slidr_cluster_tails.py"
    )
    slc_cols = extract_tsv_header_from_file(slc_path, r'out\.write\("([^\n]+?\\t[^\n]+?)\\n"\)')
    expected_slc = contracts["cluster_summary.tsv"]["columns"]
    assert slc_cols == expected_slc, (
        f"Schema drift in {slc_path}: code columns {slc_cols} != contract columns {expected_slc}"
    )
    print("✓ SLIDR_CLUSTERTAILS (cluster_summary.tsv) columns match live code.")

    # 7. SLIDR_CONSENSUS (sl_summary.tsv and sl_site_pool.tsv)
    slcon_path = (
        REPO_ROOT
        / "modules"
        / "local"
        / "slidr"
        / "consensus"
        / "resources"
        / "usr"
        / "bin"
        / "slidr_consensus.py"
    )
    slcon_text = slcon_path.read_text(encoding="utf-8")
    m_pool = re.search(r'out\.write\("([^\n]+?\\t[^\n]+?)"\s*\n\s*"([^\n]+?\\t[^\n]+?)\\n"\)', slcon_text)
    assert m_pool, f"Could not extract pool header from {slcon_path}"
    pool_cols = [c.strip() for c in (m_pool.group(1) + m_pool.group(2)).split(r"\t") if c.strip()]
    expected_pool = contracts["sl_site_pool.tsv"]["columns"]
    assert pool_cols == expected_pool, (
        f"Schema drift in {slcon_path} (pool): {pool_cols} != {expected_pool}"
    )

    m_sum = re.search(r'out\.write\("([^\n]+?\\t[^\n]+?)"\s*\n\s*"([^\n]+?\\t[^\n]+?)\\n"\)', slcon_text[m_pool.end():])
    assert m_sum, f"Could not extract summary header from {slcon_path}"
    sum_cols = [c.strip() for c in (m_sum.group(1) + m_sum.group(2)).split(r"\t") if c.strip()]
    expected_sum = contracts["sl_summary.tsv"]["columns"]
    assert sum_cols == expected_sum, (
        f"Schema drift in {slcon_path} (summary): {sum_cols} != {expected_sum}"
    )
    print("✓ SLIDR_CONSENSUS (sl_summary.tsv & sl_site_pool.tsv) columns match live code.")

    # 8. SLIDR_EXTRACTTAILS (acceptor_sites.tsv.gz and tail_stats.json)
    ext_path = (
        REPO_ROOT
        / "modules"
        / "local"
        / "slidr"
        / "extracttails"
        / "resources"
        / "usr"
        / "bin"
        / "slidr_extract_tails.py"
    )
    ext_cols = extract_tsv_header_from_file(ext_path, r'tsv_out\.write\("([^\n]+?\\t[^\n]+?)\\n"\)')
    expected_ext = contracts["acceptor_sites.tsv.gz"]["columns"]
    assert ext_cols == expected_ext, (
        f"Schema drift in {ext_path}: {ext_cols} != {expected_ext}"
    )
    ext_content = ext_path.read_text(encoding="utf-8")
    for req_key in contracts["tail_stats.json"]["required_keys"]:
        assert f'"{req_key}"' in ext_content, f"Missing stats key '{req_key}' in {ext_path}"
    print("✓ SLIDR_EXTRACTTAILS (acceptor_sites.tsv.gz & tail_stats.json) match live code.")

    # 9. SLIDR_SLRNAFILTER (sl_rna_candidates.tsv)
    slrn_path = (
        REPO_ROOT
        / "modules"
        / "local"
        / "slidr"
        / "slrnafilter"
        / "resources"
        / "usr"
        / "bin"
        / "slidr_filter_sl_rna.py"
    )
    slrn_cols = extract_tsv_header_from_file(slrn_path, r'tsv_out\.write\("([^\n]+?\\t[^\n]+?)\\n"\)')
    expected_slrn = contracts["sl_rna_candidates.tsv"]["columns"]
    assert slrn_cols == expected_slrn, (
        f"Schema drift in {slrn_path}: {slrn_cols} != {expected_slrn}"
    )
    print("✓ SLIDR_SLRNAFILTER (sl_rna_candidates.tsv) columns match live code.")

    # 10. SLIDR_TSSEXPORTER (tss_metrics.tsv)
    tss_path = (
        REPO_ROOT
        / "modules"
        / "local"
        / "slidr"
        / "tssexporter"
        / "resources"
        / "usr"
        / "bin"
        / "slidr_tss_exporter.py"
    )
    tss_cols = extract_tsv_header_from_file(tss_path, r'ftsv\.write\("([^\n]+?\\t[^\n]+?)\\n"\)')
    expected_tss = contracts["tss_metrics.tsv"]["columns"]
    assert tss_cols == expected_tss, (
        f"Schema drift in {tss_path}: {tss_cols} != {expected_tss}"
    )
    print("✓ SLIDR_TSSEXPORTER (tss_metrics.tsv) columns match live code.")


def parse_stub_cols(raw: str) -> list[str]:
    cleaned = re.sub(r"\\+n$", "", raw).strip()
    return [c.strip() for c in re.split(r"\\+t|\t", cleaned) if c.strip()]


def verify_module_stub_contracts(tsv_contracts: dict) -> None:
    """Verify stub blocks in modules emit schema-compliant headers and files."""
    contracts = tsv_contracts["contracts"]

    # 1. SELECTBESTMODEL
    sbm_main = REPO_ROOT / "modules" / "local" / "selectbestmodel" / "main.nf"
    sbm_stub = sbm_main.read_text(encoding="utf-8").split("stub:")[-1]
    m_sbm = re.search(r'printf\s+"([^"]+)"\s*>\s*"\$\{prefix\}_decision\.tsv"', sbm_stub)
    assert m_sbm, f"Missing printf stub header in {sbm_main}"
    sbm_stub_cols = parse_stub_cols(m_sbm.group(1))
    assert sbm_stub_cols == contracts["tournament_decision.tsv"]["columns"], (
        f"Stub mismatch in {sbm_main}: {sbm_stub_cols} != {contracts['tournament_decision.tsv']['columns']}"
    )
    m_sbm_busco = re.search(r'printf\s+"([^"]+)"\s*>\s*"\$\{prefix\}_busco_full_table\.tsv"', sbm_stub)
    assert m_sbm_busco, f"Missing printf stub header for busco full table in {sbm_main}"
    sbm_busco_stub_cols = parse_stub_cols(m_sbm_busco.group(1))
    assert sbm_busco_stub_cols == contracts["busco_full_table.tsv"]["columns"], (
        f"Stub mismatch in {sbm_main}: {sbm_busco_stub_cols} != {contracts['busco_full_table.tsv']['columns']}"
    )

    # 2. TESCORER
    tesc_main = REPO_ROOT / "modules" / "local" / "tescorer" / "main.nf"
    tesc_stub = tesc_main.read_text(encoding="utf-8").split("stub:")[-1]
    m_tesc = re.search(r'printf\s+"([^"]+)"\s*>\s*"\$\{prefix\}\.gene_te_scores\.tsv"', tesc_stub)
    assert m_tesc, f"Missing printf stub header in {tesc_main}"
    tesc_stub_cols = parse_stub_cols(m_tesc.group(1))
    assert tesc_stub_cols == contracts["gene_te_scores.tsv"]["columns"], (
        f"Stub mismatch in {tesc_main}: {tesc_stub_cols} != {contracts['gene_te_scores.tsv']['columns']}"
    )

    # 3. SATELLITEFILTER
    satf_main = REPO_ROOT / "modules" / "local" / "satellitefilter" / "main.nf"
    satf_stub = satf_main.read_text(encoding="utf-8").split("stub:")[-1]
    m_satf = re.search(r'printf\s+"([^"]+)"\s*>\s*\$\{prefix\}_trf_summary\.tsv', satf_stub)
    assert m_satf, f"Missing printf stub header in {satf_main}"
    satf_stub_cols = parse_stub_cols(m_satf.group(1))
    assert satf_stub_cols == contracts["trf_summary.tsv"]["columns"], (
        f"Stub mismatch in {satf_main}: {satf_stub_cols} != {contracts['trf_summary.tsv']['columns']}"
    )

    # 4. SATELLITEPROFILER
    satp_main = REPO_ROOT / "modules" / "local" / "satelliteprofiler" / "main.nf"
    satp_stub = satp_main.read_text(encoding="utf-8").split("stub:")[-1]
    m_satp = re.search(r'printf\s+"([^"]+)"\s*>\s*\$\{prefix\}_satellite_families\.tsv', satp_stub)
    assert m_satp, f"Missing printf stub header in {satp_main}"
    satp_stub_cols = parse_stub_cols(m_satp.group(1))
    assert satp_stub_cols == contracts["satellite_families.tsv"]["columns"], (
        f"Stub mismatch in {satp_main}: {satp_stub_cols} != {contracts['satellite_families.tsv']['columns']}"
    )

    # 5. SLOPPR_OPERONS
    slop_main = REPO_ROOT / "modules" / "local" / "sloppr" / "operons" / "main.nf"
    slop_stub = slop_main.read_text(encoding="utf-8").split("stub:")[-1]
    m_slop = re.search(r'printf\s+"([^"]+)"\s*>\s*\$\{prefix\}_operons_summary\.tsv', slop_stub)
    assert m_slop, f"Missing printf stub header in {slop_main}"
    slop_stub_cols = parse_stub_cols(m_slop.group(1))
    assert slop_stub_cols == contracts["operons_summary.tsv"]["columns"], (
        f"Stub mismatch in {slop_main}: {slop_stub_cols} != {contracts['operons_summary.tsv']['columns']}"
    )

    # 6. SLIDR_CLUSTERTAILS
    slc_main = REPO_ROOT / "modules" / "local" / "slidr" / "clustertails" / "main.nf"
    slc_stub = slc_main.read_text(encoding="utf-8").split("stub:")[-1]
    m_slc = re.search(r'printf\s+[\'"]([^\'"]+)[\'"]\s*>\s*\$\{prefix\}_cluster_summary\.tsv', slc_stub)
    assert m_slc, f"Missing printf stub header in {slc_main}"
    slc_stub_cols = parse_stub_cols(m_slc.group(1))
    assert slc_stub_cols == contracts["cluster_summary.tsv"]["columns"], (
        f"Stub mismatch in {slc_main}: {slc_stub_cols} != {contracts['cluster_summary.tsv']['columns']}"
    )

    # 7. SLIDR_CONSENSUS
    slcon_main = REPO_ROOT / "modules" / "local" / "slidr" / "consensus" / "main.nf"
    slcon_stub = slcon_main.read_text(encoding="utf-8").split("stub:")[-1]
    m_slsum = re.search(r'printf\s+"([^"]+)"\s*>\s*\$\{prefix\}_sl_summary\.tsv', slcon_stub)
    assert m_slsum, f"Missing printf stub header for sl_summary in {slcon_main}"
    slsum_stub_cols = parse_stub_cols(m_slsum.group(1))
    assert slsum_stub_cols == contracts["sl_summary.tsv"]["columns"], (
        f"Stub mismatch in {slcon_main}: {slsum_stub_cols} != {contracts['sl_summary.tsv']['columns']}"
    )

    # 8. SLIDR_EXTRACTTAILS
    ext_main = REPO_ROOT / "modules" / "local" / "slidr" / "extracttails" / "main.nf"
    ext_stub = ext_main.read_text(encoding="utf-8").split("stub:")[-1]
    m_ext = re.search(r'printf\s+"([^"]+)"\s*\|\s*gzip\s+-c\s*>\s*\$\{prefix\}\.acceptor_sites\.tsv\.gz', ext_stub)
    assert m_ext, f"Missing gzip stub header for acceptor_sites in {ext_main}"
    ext_stub_cols = parse_stub_cols(m_ext.group(1))
    assert ext_stub_cols == contracts["acceptor_sites.tsv.gz"]["columns"], (
        f"Stub mismatch in {ext_main}: {ext_stub_cols} != {contracts['acceptor_sites.tsv.gz']['columns']}"
    )

    # 9. SLIDR_SLRNAFILTER
    slrn_main = REPO_ROOT / "modules" / "local" / "slidr" / "slrnafilter" / "main.nf"
    slrn_stub = slrn_main.read_text(encoding="utf-8").split("stub:")[-1]
    m_slrn = re.search(r'printf\s+"([^"]+)"\s*>\s*\$\{prefix\}\.sl_rna_candidates\.tsv', slrn_stub)
    assert m_slrn, f"Missing printf stub header in {slrn_main}"
    slrn_stub_cols = parse_stub_cols(m_slrn.group(1))
    assert slrn_stub_cols == contracts["sl_rna_candidates.tsv"]["columns"], (
        f"Stub mismatch in {slrn_main}: {slrn_stub_cols} != {contracts['sl_rna_candidates.tsv']['columns']}"
    )

    # 10. SLIDR_TSSEXPORTER
    tss_main = REPO_ROOT / "modules" / "local" / "slidr" / "tssexporter" / "main.nf"
    tss_stub = tss_main.read_text(encoding="utf-8").split("stub:")[-1]
    m_tss = re.search(r'printf\s+"([^"]+)"\s*>\s*\$\{prefix\}\.tss_metrics\.tsv', tss_stub)
    assert m_tss, f"Missing printf stub header in {tss_main}"
    tss_stub_cols = parse_stub_cols(m_tss.group(1))
    assert tss_stub_cols == contracts["tss_metrics.tsv"]["columns"], (
        f"Stub mismatch in {tss_main}: {tss_stub_cols} != {contracts['tss_metrics.tsv']['columns']}"
    )

    print("✓ All module stub contracts emit valid schema-compliant headers and deliverables.")


def verify_extract_and_stamp_tools() -> None:
    """Verify that extract_code_contracts.py --check and stamp_module_shas.py --check exit 0."""
    extractor = REPO_ROOT / "scripts" / "extract_code_contracts.py"
    stamper = REPO_ROOT / "scripts" / "stamp_module_shas.py"

    res_ext = subprocess.run([sys.executable, str(extractor), "--check"], cwd=REPO_ROOT, capture_output=True, text=True)
    assert res_ext.returncode == 0, f"extract_code_contracts.py --check failed:\n{res_ext.stderr}\n{res_ext.stdout}"

    res_st = subprocess.run([sys.executable, str(stamper), "--check"], cwd=REPO_ROOT, capture_output=True, text=True)
    assert res_st.returncode == 0, f"stamp_module_shas.py --check failed:\n{res_st.stderr}\n{res_st.stdout}"

    print("✓ extract_code_contracts.py --check and stamp_module_shas.py --check passed with exit code 0.")


def main() -> int:
    tsv_contracts = verify_schemas_structure()
    verify_live_code_contracts(tsv_contracts)
    verify_module_stub_contracts(tsv_contracts)
    verify_extract_and_stamp_tools()
    print("check_data_contracts: PASS (all TSV/JSON schemas, live code contracts, and stubs verified)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
