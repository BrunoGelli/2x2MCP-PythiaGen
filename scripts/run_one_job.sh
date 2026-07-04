#!/usr/bin/env bash
set -euo pipefail
JOB_LIST=${JOB_LIST:-job_list.csv}; JOB_ID=${1:-${SLURM_ARRAY_TASK_ID:-}}; EXEC=${EXEC:-./allInOne_multimeson}; NTHREADS=${NTHREADS:-1}; OVERWRITE=${OVERWRITE:-0}
[[ -n "$JOB_ID" ]] || { echo "Usage: $0 JOB_ID"; exit 2; }
row=$(python3 - "$JOB_LIST" "$JOB_ID" <<'PY'
import csv,sys
for r in csv.DictReader(open(sys.argv[1])):
    if r['job_id']==str(sys.argv[2]):
        print('|'.join(r[k] for k in ['seed','n_events','mcp_mass_GeV','geometry','emitter_name','production_mode','output_file','job_id'])); break
else: sys.exit(1)
PY
)
IFS='|' read -r seed nevents mass geom emitter mode outfile jid <<< "$row"
mkdir -p "$(dirname "$outfile")" logs/jobs
if [[ -s "$outfile" && "$OVERWRITE" != 1 ]]; then echo "Output exists, skipping: $outfile"; exit 0; fi
log="logs/jobs/job_${jid}.log"
cmd=("$EXEC" "$seed" "$NTHREADS" "$nevents" "$mass" "$geom" "$emitter" "$mode" "$outfile" --mode fixed-events --n-events "$nevents" --write-spectra "${WRITE_SPECTRA:-accepted}" --spectra-prescale "${SPECTRA_PRESCALE:-1}" --batch)
if [[ -n "${PRODUCTION_CONFIG:-}" ]]; then cmd+=(--production-config "$PRODUCTION_CONFIG"); fi
"${cmd[@]}" >"$log" 2>&1
