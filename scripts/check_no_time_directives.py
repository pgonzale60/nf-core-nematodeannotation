#!/usr/bin/env python3
"""Pre-commit check to enforce that Nextflow configuration files never specify 'time' directives.

On Farm22 LSF clusters, walltime is governed strictly by the LSF queue itself
(e.g., normal <= 12h, long <= 48h, basement). Declaring explicit `time = ...`
in process configuration directives triggers artificial LSF RUNLIMIT rejection errors.
"""

import pathlib
import re
import sys

REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent

# Matches `time = ...` or `time: ...` in configuration lines
TIME_DIRECTIVE_RE = re.compile(r"^\s*time\s*[:=]\s*.*$", re.MULTILINE)


def check_config_file(filepath: pathlib.Path) -> list[tuple[int, str]]:
    violations = []
    text = filepath.read_text(encoding="utf-8")
    for idx, line in enumerate(text.splitlines(), start=1):
        # Ignore comments
        stripped = line.strip()
        if stripped.startswith("//") or stripped.startswith("/*") or stripped.startswith("*"):
            continue
        if re.match(r"^\s*time\s*[:=]", line):
            violations.append((idx, line.strip()))
    return violations


def main() -> int:
    config_files = list(REPO_ROOT.glob("conf/**/*.config")) + list(REPO_ROOT.glob("*.config"))
    failed = False

    for cfg in sorted(config_files):
        violations = check_config_file(cfg)
        if violations:
            failed = True
            print(f"ERROR: Found prohibited 'time' directive in {cfg.relative_to(REPO_ROOT)}:")
            for line_no, content in violations:
                print(f"  Line {line_no}: {content}")

    if failed:
        print("\nRule Violation: Nextflow configs must NEVER specify 'time' directives.")
        print("Walltime is governed strictly and solely by queue selection on Farm22 LSF.")
        print("Please remove all 'time' directives from the config.\n")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
