#!/usr/bin/env python3

# Checks that every tracked notebook is run by a Snakemake rule and has its executed copy

from __future__ import annotations

import os
import re
import subprocess
import sys
from pathlib import Path

import yaml

REPO = Path(__file__).resolve().parents[2]
SNAKEFILE = REPO / "workflow" / "Snakefile"
EXECUTED = ".executed.ipynb"
ALLOWED_EXTRA = {".gitkeep"}
ASSIGNMENT = re.compile(r"^([A-Z][A-Z0-9_]*)\s*=\s*(.+)$")
NOTEBOOK_LITERAL = re.compile(r"""(f?)(["'])([^"'\n]*\.ipynb)\2""")
PLACEHOLDER = re.compile(r"\{([A-Za-z_][A-Za-z0-9_]*)\}")


# Merges one config mapping into another, recursively.
def merge(base: dict, extra: dict) -> dict:
    for key, value in extra.items():
        if isinstance(value, dict) and isinstance(base.get(key), dict):
            merge(base[key], value)
        else:
            base[key] = value
    return base


# Lists the config files and rule files the Snakefile loads, in order.
def workflow_files() -> tuple[list[Path], list[Path]]:
    text = SNAKEFILE.read_text()
    configs = [REPO / p for p in re.findall(r'^configfile:\s*"([^"]+)"', text, re.M)]
    rules = [SNAKEFILE.parent / p for p in re.findall(r'^include:\s*"([^"]+)"', text, re.M)]
    return configs, rules


# Yields every string value in a nested config.
def config_strings(node):
    if isinstance(node, dict):
        for value in node.values():
            yield from config_strings(value)
    elif isinstance(node, list):
        for value in node:
            yield from config_strings(value)
    elif isinstance(node, str):
        yield node


# Resolves the source notebook paths referenced by the workflow.
def referenced_notebooks() -> tuple[set[str], list[str]]:
    configs, rules = workflow_files()
    config: dict = {}
    for path in configs:
        merge(config, yaml.safe_load(path.read_text()) or {})

    refs = {s for s in config_strings(config) if s.endswith(".ipynb") and not s.endswith(EXECUTED)}
    unresolved = []
    names = {"config": config, "os": os, "REPO_ROOT": ""}
    for path in rules:
        for line in path.read_text().splitlines():
            match = ASSIGNMENT.match(line)
            if match:
                try:
                    names[match.group(1)] = eval(match.group(2), {"__builtins__": {}}, names)
                except Exception:
                    pass
            for is_fstring, _, literal in NOTEBOOK_LITERAL.findall(line):
                if literal.endswith(EXECUTED):
                    continue
                if is_fstring:
                    literal = PLACEHOLDER.sub(
                        lambda m: str(names.get(m.group(1), m.group(0))), literal
                    )
                if "{" in literal or "/" not in literal:
                    unresolved.append(f"{path.relative_to(REPO)}: {literal}")
                    continue
                refs.add(os.path.normpath(literal))
    return refs, unresolved


# Lists files tracked by git under nbs/.
def tracked_files() -> list[str]:
    out = subprocess.run(
        ["git", "ls-files", "nbs"], cwd=REPO, check=True, capture_output=True, text=True
    ).stdout
    return [line for line in out.splitlines() if line]


def main() -> int:
    tracked = tracked_files()
    tracked_set = set(tracked)
    sources = {p for p in tracked if p.endswith(".ipynb") and not p.endswith(EXECUTED)}
    executed = {p for p in tracked if p.endswith(EXECUTED)}
    refs, unresolved = referenced_notebooks()
    errors = [f"cannot resolve notebook path in {item}" for item in unresolved]

    for path in sorted(refs):
        if not path.startswith("nbs/"):
            continue
        if path not in tracked_set:
            errors.append(f"{path}: referenced by the workflow but not tracked")
        elif path.removesuffix(".ipynb") + EXECUTED not in tracked_set:
            errors.append(f"{path}: run by Snakemake but its {EXECUTED} copy is not tracked")

    for path in sorted(sources - refs):
        errors.append(f"{path}: not run by any Snakemake rule")

    for path in sorted(executed):
        if path.removesuffix(EXECUTED) + ".ipynb" not in tracked_set:
            errors.append(f"{path}: executed copy without its source notebook")

    for path in tracked:
        name = os.path.basename(path)
        if not name.endswith(".ipynb") and name not in ALLOWED_EXTRA:
            errors.append(f"{path}: only notebooks belong in nbs/")

    if errors:
        print(f"{len(errors)} notebook problem(s):", file=sys.stderr)
        for error in errors:
            print(f"  {error}", file=sys.stderr)
        return 1
    print(f"OK: {len(sources)} notebooks, all run by Snakemake")
    return 0


if __name__ == "__main__":
    sys.exit(main())
