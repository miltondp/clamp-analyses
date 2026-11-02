#!/usr/bin/env python3
"""Check the execution environment without accessing a workspace or running GLS."""
import argparse
from importlib.metadata import distribution, version
import json
from pathlib import Path
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--require-version", default="0.7.3")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    # Import the packages used by collection too: distribution metadata alone
    # cannot detect a broken numerical-library installation.
    import phenoplier
    import numpy
    import pandas
    import scipy
    import h5py

    if phenoplier.__version__ != args.require_version:
        parser.error(f"phenoplier {phenoplier.__version__} installed; require {args.require_version}")
    source = json.loads(distribution("phenoplier").read_text("direct_url.json") or "{}")
    report = dict(python=sys.version.split()[0], executable=sys.executable,
                  phenoplier=phenoplier.__version__,
                  phenoplier_commit=source.get("vcs_info", {}).get("commit_id"),
                  numpy=numpy.__version__,
                  pandas=pandas.__version__, scipy=scipy.__version__, h5py=h5py.__version__,
                  snakemake=version("snakemake"), statsmodels=version("statsmodels"))
    text = json.dumps(report, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()
