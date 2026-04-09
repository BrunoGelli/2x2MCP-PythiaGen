# Multi-meson MCP production + acceptance pipeline (updated)

## Main changes in this revision

This revision turns the generator into a more production-friendly tool:

- Pythia's own periodic `Pythia::next()` chatter is explicitly disabled.
- The accepted-event target is now **global per mass point**, not per thread.
- A **wall-time safety limit** can stop a mass point before it runs forever.
- The progress monitor is more informative:
  - accepted progress bar
  - generated-event count
  - `pythia.next()` failure count
  - event rate
  - accepted-event rate
  - elapsed time
  - ETA
  - stop reason
- The ROOT summary tree now stores:
  - `n_next_calls`
  - `n_next_failures`
  - `stop_reason`
  - `stop_reason_name`
- The code has been more heavily commented around the physics-sensitive logic.

## Important behavior change

The third command-line argument is now:

```text
acceptedTargetTotal
```

This is the **total accepted-event target for the whole mass point**, shared across all worker threads.

So if you run with 8 threads and `acceptedTargetTotal = 500`, the job stops when the mass point has accumulated about 500 accepted MCPs in total, not 500 per thread.

This makes the generator much more natural for future production modes, especially when you later parallelize over masses at the batch level.

## New optional wall-time limit

A new optional final argument was added:

```text
[maxWallSeconds=0]
```

- `0` means disabled
- any positive value gives a hard wall-time safety stop for the whole mass point

Example:

```bash
./allInOne_multimeson 1000 8 500 0.02 2x2 all 1000000 mcp_acceptance_results.root beam.config momentum.config 1800
```

That would stop after 30 minutes if the accepted target was not reached first.

## Compile

```bash
g++ -O2 -std=c++17 -pthread allInOne_multimeson.cc -o allInOne_multimeson \
  $(root-config --cflags --libs) \
  $(pythia8-config --cxxflags --libs)
```

## Run

Minimal example:

```bash
./allInOne_multimeson 1000 8 500 0.02 2x2
```

Full example with wall-time safety:

```bash
./allInOne_multimeson 1000 8 500 0.02 2x2 all 1000000 mcp_acceptance_results.root beam.config momentum.config 1800
```

Arguments:

```text
[seed] [nThreads] [acceptedTargetTotal] [mcpMass_GeV] [geometry]
[parentList=all] [maxEventsPerThread=1000000] [outFile=mcp_acceptance_results.root]
[beamConfig=beam.config] [momentumConfig=momentum.config] [maxWallSeconds=0]
```

## Batch scan

```bash
./run_all_masses_multimeson.sh
```

Environment overrides:

```bash
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
```

## Output tree schema

Tree name: `mcp_summary`

Branches:

- `thread_id`
- `seed`
- `geometry_id`
- `geometry_name`
- `parent_pdg`
- `parent_type`
- `parent_name`
- `parent_type_name`
- `stop_reason`
- `stop_reason_name`
- `mcp_mass`
- `n_events_generated`
- `n_next_calls`
- `n_next_failures`
- `n_parent_total`
- `n_mcp_total`
- `n_mcp_accepted`
- `acceptance_fraction`
- `parent_yield_per_event`

Each row is still one `(thread, mass, parent)` summary row.

## Plotting

The existing plotting macro should continue to work because the original branches it uses are still present.

Inside ROOT:

```cpp
.x plot_acceptance_vs_mass_multimeson.C
plot_acceptance_vs_mass_multimeson("mcp_acceptance_results.root", 1.5e19, 1.0e-2);
```

## One important unit reminder

The command-line mass is expected in **GeV**.

Examples:

- `10 MeV -> 0.01 GeV`
- `1000 MeV -> 1.0 GeV`
