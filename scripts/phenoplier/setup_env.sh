#!/usr/bin/env bash
# Install the pinned PhenoPLIER environment.
set -euo pipefail

REF="${PHENOPLIER_REF:-v0.5.2}"
CLONE_DIR="${PHENOPLIER_CLONE_DIR:-$HOME/phenoplier-cli}"
TARGET_ENV="${PHENOPLIER_TARGET_ENV:-phenoplier-cli-neo}"

curl -fsSL "https://raw.githubusercontent.com/pivlab/phenoplier-cli/${REF}/scripts/install_phenoplier_cli.sh" \
  | bash -s -- --dir "$CLONE_DIR" --env "$TARGET_ENV" --ref "$REF" --extras slurm "$@"

cat <<EOF

Next: wire the reference data bundle into the workspace once:
  conda activate $TARGET_ENV
  phenoplier workspace init
  phenoplier workspace link /path/to/phenoplier_full_data
EOF
