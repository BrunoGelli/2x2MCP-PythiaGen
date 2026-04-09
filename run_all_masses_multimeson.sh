#!/bin/bash
set -euo pipefail

NTHREADS="${NTHREADS:-8}"
NACCEPT_TOTAL="${NACCEPT_TOTAL:-500}"
MAX_EVENTS_PER_THREAD="${MAX_EVENTS_PER_THREAD:-1000000}"
MAX_WALL_SECONDS="${MAX_WALL_SECONDS:-0}"
SEED_BASE="${SEED_BASE:-1000}"
EXEC="${EXEC:-./allInOne_multimeson}"
ROOT_FILE="${ROOT_FILE:-mcp_acceptance_results.root}"
GEOMETRY="${GEOMETRY:-2x2}"
PARENTS="${PARENTS:-all}"
BEAM_CONFIG="${BEAM_CONFIG:-beam.config}"
MOMENTUM_CONFIG="${MOMENTUM_CONFIG:-momentum.config}"
MASSES_FILE="${MASSES_FILE:-masses.txt}"

usage() {
  cat <<USAGE
Usage:
  $0 [masses_file]

Environment overrides:
  NTHREADS=8
  NACCEPT_TOTAL=500
  MAX_EVENTS_PER_THREAD=1000000
  MAX_WALL_SECONDS=0
  SEED_BASE=1000
  EXEC=./allInOne_multimeson
  ROOT_FILE=mcp_acceptance_results.root
  GEOMETRY=2x2
  PARENTS=all
  BEAM_CONFIG=beam.config
  MOMENTUM_CONFIG=momentum.config
  MASSES_FILE=masses.txt

Notes:
  - NACCEPT_TOTAL is the total accepted-event target per mass point.
  - MAX_WALL_SECONDS=0 disables the wall-time safety stop.
USAGE
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  usage
  exit 0
fi

if [[ $# -ge 1 ]]; then
  MASSES_FILE="$1"
fi

if [[ ! -f "$MASSES_FILE" ]]; then
  echo "Error: masses file '$MASSES_FILE' not found."
  exit 1
fi

if [[ ! -x "$EXEC" ]]; then
  echo "Error: executable '$EXEC' not found or not executable."
  exit 1
fi

mapfile -t RAW_MASSES < "$MASSES_FILE"
MASSES=()
for M in "${RAW_MASSES[@]}"; do
  M="$(echo "$M" | xargs)"
  [[ -z "$M" ]] && continue
  [[ "$M" =~ ^# ]] && continue
  MASSES+=("$M")
done

if [[ ${#MASSES[@]} -eq 0 ]]; then
  echo "Error: no masses found in '$MASSES_FILE'."
  exit 1
fi

rm -f "$ROOT_FILE"

draw_overall_bar() {
  local current="$1"
  local total="$2"
  local width=36
  local filled=$(( current * width / total ))
  local empty=$(( width - filled ))
  local bar
  printf -v bar '%*s' "$filled" ''
  bar=${bar// /#}
  local rest
  printf -v rest '%*s' "$empty" ''
  rest=${rest// /-}
  printf '[%s%s] %d/%d' "$bar" "$rest" "$current" "$total"
}

echo "Running MCP scan"
echo "  executable:             $EXEC"
echo "  geometry:               $GEOMETRY"
echo "  parents:                $PARENTS"
echo "  threads per mass:       $NTHREADS"
echo "  accepted target / mass: $NACCEPT_TOTAL"
echo "  max events / thread:    $MAX_EVENTS_PER_THREAD"
echo "  max wall time [s]:      $MAX_WALL_SECONDS"
echo "  output root file:       $ROOT_FILE"
echo "  masses file:            $MASSES_FILE"
echo

TOTAL_MASSES=${#MASSES[@]}
START_TS=$(date +%s)

for (( MASS_INDEX=0; MASS_INDEX<TOTAL_MASSES; ++MASS_INDEX )); do
  MASS="${MASSES[$MASS_INDEX]}"
  SEED=$((SEED_BASE + MASS_INDEX * 1000))

  echo "======================================================================"
  echo "Scan progress: $(draw_overall_bar "$MASS_INDEX" "$TOTAL_MASSES")"
  echo "Starting mass $(printf '%s' "$MASS") GeV   |   seed $SEED"

  "$EXEC" \
    "$SEED" \
    "$NTHREADS" \
    "$NACCEPT_TOTAL" \
    "$MASS" \
    "$GEOMETRY" \
    "$PARENTS" \
    "$MAX_EVENTS_PER_THREAD" \
    "$ROOT_FILE" \
    "$BEAM_CONFIG" \
    "$MOMENTUM_CONFIG" \
    "$MAX_WALL_SECONDS"

done

END_TS=$(date +%s)
ELAPSED=$((END_TS - START_TS))

echo "======================================================================"
echo "Scan progress: $(draw_overall_bar "$TOTAL_MASSES" "$TOTAL_MASSES")"
echo "Done. Results written to $ROOT_FILE"
echo "Elapsed: ${ELAPSED}s"
