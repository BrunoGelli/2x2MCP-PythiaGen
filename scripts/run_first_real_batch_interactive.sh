#!/bin/bash
set -euo pipefail

EXEC="${EXEC:-./allInOne_multimeson}"
MASSES_FILE="${MASSES_FILE:-masses_meson_first_real_batch.txt}"
RESULTS_DIR="${RESULTS_DIR:-batch_outputs_$(date +%Y%m%d_%H%M%S)}"
GEOMETRY="${GEOMETRY:-2x2}"
PARENTS="${PARENTS:-all}"
BEAM_CONFIG="${BEAM_CONFIG:-beam.config}"
MOMENTUM_CONFIG="${MOMENTUM_CONFIG:-momentum.config}"
SEED_BASE="${SEED_BASE:-50000}"
MAX_PARALLEL="${MAX_PARALLEL:-8}"
ACCEPT_TARGET_TOTAL="${ACCEPT_TARGET_TOTAL:-300}"
WALL_SECONDS="${WALL_SECONDS:-1800}"

if [[ ! -x "$EXEC" ]]; then
  echo "Error: executable '$EXEC' not found or not executable."
  return 1 2>/dev/null || exit 1
fi

if [[ ! -f "$MASSES_FILE" ]]; then
  echo "Error: masses file '$MASSES_FILE' not found."
  return 1 2>/dev/null || exit 1
fi

mkdir -p "$RESULTS_DIR"/logs

mapfile -t RAW_MASSES < "$MASSES_FILE"
MASSES=()
for M in "${RAW_MASSES[@]}"; do
  M="$(echo "$M" | xargs)"
  [[ -z "$M" ]] && continue
  [[ "$M" =~ ^# ]] && continue
  MASSES+=("$M")
done

choose_threads() {
  local mass="$1"
  python - "$mass" <<'PY'
import sys
m = float(sys.argv[1])
if m < 0.10:
    print(2)
elif m < 0.40:
    print(4)
else:
    print(12)
PY
}

choose_max_events_per_thread() {
  local mass="$1"
  python - "$mass" <<'PY'
import sys
m = float(sys.argv[1])
if m < 0.10:
    print(250000)
elif m < 0.40:
    print(1000000)
else:
    print(5000000)
PY
}

launch_one() {
  local idx="$1"
  local mass="$2"
  local seed=$((SEED_BASE + idx * 1000))
  local nthreads
  local max_events
  nthreads="$(choose_threads "$mass")"
  max_events="$(choose_max_events_per_thread "$mass")"

  local mass_tag
  mass_tag="$(printf "%0.6f" "$mass")"
  local root_out="${RESULTS_DIR}/mcp_mass_${mass_tag}.root"
  local log_out="${RESULTS_DIR}/logs/mcp_mass_${mass_tag}.log"

  echo "Launching mass=${mass_tag} GeV  threads=${nthreads}  maxEvents/thread=${max_events}  seed=${seed}"

  "$EXEC"     "$seed"     "$nthreads"     "$ACCEPT_TARGET_TOTAL"     "$mass"     "$GEOMETRY"     "$PARENTS"     "$max_events"     "$root_out"     "$BEAM_CONFIG"     "$MOMENTUM_CONFIG"     "$WALL_SECONDS"     > "$log_out" 2>&1
}

echo "Interactive batch configuration"
echo "  executable:              $EXEC"
echo "  masses file:             $MASSES_FILE"
echo "  results dir:             $RESULTS_DIR"
echo "  max parallel processes:  $MAX_PARALLEL"
echo "  accepted target / mass:  $ACCEPT_TARGET_TOTAL"
echo "  wall time / mass [s]:    $WALL_SECONDS"
echo

active=0
for i in "${!MASSES[@]}"; do
  launch_one "$i" "${MASSES[$i]}" &
  active=$((active + 1))
  if [[ "$active" -ge "$MAX_PARALLEL" ]]; then
    wait -n
    active=$((active - 1))
  fi
done

wait

echo
echo "All mass points finished."
echo "Outputs are in: $RESULTS_DIR"
echo "Merge later with:"
echo "  hadd -f ${RESULTS_DIR}/mcp_acceptance_merged.root ${RESULTS_DIR}/mcp_mass_*.root"
