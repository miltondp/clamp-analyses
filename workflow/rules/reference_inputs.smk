# Pinned small references: shipped snapshots or checksum-verified downloads.
import json
import re

with open("workflow/config/reference_inputs.json") as handle:
    REFERENCE_ASSETS = json.load(handle)["assets"]
REFERENCE_PATHS = {asset["path"]: key for key, asset in REFERENCE_ASSETS.items()}
REFERENCE_DOWNLOAD_PATTERN = "|".join(
    re.escape(asset["path"]) for asset in REFERENCE_ASSETS.values()
    if not asset.get("existing_rule")
)

rule fetch_reference_input:
    output:
        "{reference_file}",
    wildcard_constraints:
        reference_file=REFERENCE_DOWNLOAD_PATTERN,
    params:
        key=lambda wc: REFERENCE_PATHS[wc.reference_file],
    conda: "snakemake"
    shell:
        "python scripts/inputs/reference_inputs.py --only {params.key:q}"

rule reference_inputs:
    input:
        [asset["path"] for asset in REFERENCE_ASSETS.values()],
