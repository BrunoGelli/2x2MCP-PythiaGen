#!/bin/bash
set -euo pipefail

RESULTS_DIR="${1:-}"
OUTFILE="${2:-}"

if [[ -z "$RESULTS_DIR" ]]; then
  echo "Usage: $0 <results_dir> [outfile]"
  exit 1
fi

if [[ ! -d "$RESULTS_DIR" ]]; then
  echo "Error: directory '$RESULTS_DIR' not found."
  exit 1
fi

if [[ -z "$OUTFILE" ]]; then
  OUTFILE="${RESULTS_DIR}/mcp_acceptance_merged.root"
fi

hadd -f "$OUTFILE" "${RESULTS_DIR}"/mcp_mass_*.root
echo "Merged output written to $OUTFILE"
