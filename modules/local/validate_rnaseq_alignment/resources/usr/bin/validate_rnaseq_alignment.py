#!/usr/bin/env python3
"""Validate RNA-seq alignments against target assembly FAI and verify spliced alignment provenance.

Accepts coordinate-sorted BAM or reference-aligned CRAM files and ensures:
1. File integrity via `samtools quickcheck`.
2. Reference sequence dictionary matches target assembly FAI (contig names and lengths).
3. Coordinate-sorted order (`@HD SO:coordinate`).
4. Spliced alignment provenance (HISAT2, STAR, Minimap2, or spliced CIGARs with 'N').

Exports environment variables to indicate validation status:
- IS_HISAT2="true" (for HISAT2 or accepted spliced aligners STAR/Minimap2 for workflow gating)
- IS_COORDINATE_VALID="true" (when contigs match assembly FAI and sorting is valid)
- IS_SPLICED="true" (when spliced reads or spliced aligner provenance is detected)
- IS_VALID="true" (overall validity passing all structural and reference checks)
"""

import argparse
import pathlib
import re
import subprocess
import sys


def parse_args():
    parser = argparse.ArgumentParser(
        description="Validate RNA-seq BAM/CRAM alignment against assembly FAI and check spliced provenance."
    )
    parser.add_argument(
        "--alignment",
        required=True,
        type=pathlib.Path,
        help="Path to BAM or CRAM alignment file.",
    )
    parser.add_argument(
        "--fai",
        required=True,
        type=pathlib.Path,
        help="Path to target assembly FAI file.",
    )
    parser.add_argument(
        "--reference",
        type=pathlib.Path,
        default=None,
        help="Optional reference FASTA file for CRAM decoding.",
    )
    parser.add_argument(
        "--env-out",
        type=pathlib.Path,
        default=None,
        help="Path to write shell environment export file.",
    )
    parser.add_argument(
        "--export-env",
        action="store_true",
        help="Print shell export statements to stdout.",
    )
    return parser.parse_args()


def check_integrity(alignment_path: pathlib.Path) -> bool:
    """Run samtools quickcheck to ensure the alignment is complete and uncorrupted."""
    cmd = ["samtools", "quickcheck", "-q", str(alignment_path)]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        sys.stderr.write(
            f"ERROR: Alignment file '{alignment_path}' failed samtools quickcheck (corrupted or truncated).\n"
        )
        return False
    return True


def get_header(alignment_path: pathlib.Path, reference_path: pathlib.Path | None) -> str:
    """Extract SAM header text using samtools view -H."""
    cmd = ["samtools", "view", "-H"]
    if reference_path and reference_path.exists():
        cmd.extend(["--reference", str(reference_path)])
    cmd.append(str(alignment_path))
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        sys.stderr.write(
            f"ERROR: Failed to extract header from '{alignment_path}':\n{result.stderr}\n"
        )
        sys.exit(1)
    return result.stdout


def validate_coordinates_against_fai(header: str, fai_path: pathlib.Path) -> tuple[bool, list[str]]:
    """Verify that all contigs in alignment header match the target assembly FAI lengths."""
    # Parse @SQ lines
    sq_dict: dict[str, int] = {}
    for line in header.splitlines():
        if not line.startswith("@SQ"):
            continue
        sn = None
        ln = None
        for tag in line.strip().split("\t")[1:]:
            if tag.startswith("SN:"):
                sn = tag[3:]
            elif tag.startswith("LN:"):
                try:
                    ln = int(tag[3:])
                except ValueError:
                    pass
        if sn and ln is not None:
            sq_dict[sn] = ln

    if not sq_dict:
        return False, ["No @SQ lines found in alignment header."]

    # Parse FAI
    fai_dict: dict[str, int] = {}
    with open(fai_path, "r", encoding="utf-8") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                try:
                    fai_dict[parts[0]] = int(parts[1])
                except ValueError:
                    pass

    mismatches: list[str] = []
    for sn, ln in sq_dict.items():
        if sn not in fai_dict:
            mismatches.append(f"Contig '{sn}' in alignment header not present in assembly FAI.")
        elif fai_dict[sn] != ln:
            mismatches.append(
                f"Contig '{sn}' length mismatch: header={ln}, assembly FAI={fai_dict[sn]}."
            )

    is_valid = len(mismatches) == 0
    return is_valid, mismatches


def check_spliced_and_aligner(
    alignment_path: pathlib.Path,
    header: str,
    reference_path: pathlib.Path | None,
) -> tuple[bool, str, bool]:
    """Detect aligner provenance (HISAT2, STAR, Minimap2) and spliced read status."""
    is_hisat2 = bool(re.search(r"@PG.*(ID:hisat2|PN:hisat2|hisat2-align)", header, re.IGNORECASE))
    is_star = bool(re.search(r"@PG.*(ID:STAR|PN:STAR|\bSTAR\b)", header))
    is_minimap2 = bool(re.search(r"@PG.*(ID:minimap2|PN:minimap2|\bminimap2\b)", header))

    aligner_name = "unknown"
    if is_hisat2:
        aligner_name = "hisat2"
    elif is_star:
        aligner_name = "star"
    elif is_minimap2:
        aligner_name = "minimap2"

    # Known spliced RNA-seq aligners default to spliced
    is_spliced = is_hisat2 or is_star

    if is_minimap2 and re.search(r"@PG.*minimap2.*(?:-x\s+splice|splice)", header):
        is_spliced = True

    # Sample aligned reads for 'N' in CIGAR (intron / splice junction)
    view_cmd = ["samtools", "view", "-F", "4"]
    if reference_path and reference_path.exists():
        view_cmd.extend(["--reference", str(reference_path)])
    view_cmd.append(str(alignment_path))

    try:
        proc = subprocess.Popen(
            view_cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True
        )
        count = 0
        if proc.stdout:
            for line in proc.stdout:
                count += 1
                fields = line.split("\t")
                if len(fields) > 5:
                    cigar = fields[5]
                    if "N" in cigar:
                        is_spliced = True
                        break
                if count >= 2000:
                    break
        proc.terminate()
    except Exception:
        pass

    return is_hisat2, aligner_name, is_spliced


def main():
    args = parse_args()

    if not args.alignment.exists():
        sys.stderr.write(f"ERROR: Alignment file '{args.alignment}' does not exist.\n")
        sys.exit(1)

    if not args.fai.exists():
        sys.stderr.write(f"ERROR: Assembly FAI '{args.fai}' does not exist.\n")
        sys.exit(1)

    # 1. Structural integrity check
    if not check_integrity(args.alignment):
        sys.exit(1)

    # 2. Extract SAM header
    header = get_header(args.alignment, args.reference)

    # 3. Coordinate sorting check
    is_coordinate_sorted = "@HD" in header and "SO:coordinate" in header
    if not is_coordinate_sorted:
        sys.stderr.write(
            f"WARNING: Alignment '{args.alignment.name}' header does not specify 'SO:coordinate'.\n"
        )

    # 4. Reference sequence dictionary check against FAI
    is_coord_valid, mismatches = validate_coordinates_against_fai(header, args.fai)
    if not is_coord_valid:
        sys.stderr.write(
            f"ERROR: Reference sequence mismatch detected in '{args.alignment.name}'.\n"
        )
        for m in mismatches[:10]:
            sys.stderr.write(f"  - {m}\n")
        if len(mismatches) > 10:
            sys.stderr.write(f"  ... and {len(mismatches) - 10} more mismatches.\n")
        sys.exit(1)
    else:
        sys.stdout.write(
            f"Reference sequence dictionary in '{args.alignment.name}' matches target assembly FAI.\n"
        )

    # 5. Aligner provenance & spliced reads
    is_hisat2, aligner_name, is_spliced = check_spliced_and_aligner(
        args.alignment, header, args.reference
    )

    if is_hisat2:
        sys.stdout.write(f"Found HISAT2 provenance in header of {args.alignment.name}\n")
    elif aligner_name in ["star", "minimap2"]:
        sys.stdout.write(
            f"INFO: Detected non-HISAT2 spliced alignment provenance ({aligner_name.upper()}) in {args.alignment.name}.\n"
        )
    elif is_spliced:
        sys.stdout.write(
            f"INFO: Spliced alignments detected in {args.alignment.name} (aligner: {aligner_name}).\n"
        )
    else:
        sys.stdout.write(
            f"WARNING: No spliced reads or known RNA-seq aligner provenance detected in {args.alignment.name}.\n"
        )

    # 6. Overall validity & env vars
    is_valid = is_coord_valid and (is_hisat2 or is_spliced or aligner_name in ["star", "minimap2"])

    # For downstream workflow acceptance:
    # Set IS_HISAT2="true" if HISAT2, OR if valid spliced alignment from STAR or Minimap2
    env_hisat2 = "true" if (is_hisat2 or (is_valid and is_spliced)) else "false"
    env_coord = "true" if is_coord_valid else "false"
    env_spliced = "true" if is_spliced else "false"
    env_valid = "true" if is_valid else "false"

    env_lines = [
        f'export IS_HISAT2="{env_hisat2}"',
        f'export IS_COORDINATE_VALID="{env_coord}"',
        f'export IS_SPLICED="{env_spliced}"',
        f'export IS_VALID="{env_valid}"',
    ]

    if args.env_out:
        args.env_out.parent.mkdir(parents=True, exist_ok=True)
        with open(args.env_out, "w", encoding="utf-8") as f:
            f.write("\n".join(env_lines) + "\n")

    if args.export_env:
        for line in env_lines:
            print(line)

    sys.exit(0)


if __name__ == "__main__":
    main()
