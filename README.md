# One-emitter MCP production + geometric acceptance pipeline

This repository estimates millicharged-particle (MCP) flux entering the DUNE ND-LAr 2x2 / LAr geometry using the factorization

`SM parent production × exotic decay × geometric acceptance × detector response`.

The generator now uses **one MCP-emitter species per run**.  For example, a `pi0` run forces only `pi0 -> gamma chi chibar`; eta, eta', omega, phi, etc. keep their Standard Model decay tables, so Pythia feed-down such as `eta -> pi0 + X` can still create secondary pi0 emitters.  Physical branching ratios and epsilon scaling are deliberately left for post-processing.

## Emitters and production modes

Central metadata lives in `mcp_emitters.h` and includes name, PDG, type, mass, default production mode, forced decay description, and placeholder branching normalization.

Supported emitters:

- pseudoscalars: `pi0` (111), `eta` (221), `etap` (331)
- vectors: `rho0` (113), `omega` (223), `phi` (333), `jpsi` (443)

Production modes:

- `light_mesons`: reads the existing beam/momentum configs and preserves the SoftQCD-style configuration already in the repo.
- `charmonium`: isolated J/psi/onium mode using separate Pythia charmonium steering.  It is recorded separately and must not be blindly normalized with light mesons.
- DY is intentionally not implemented yet; the structure reserves space for it later.

Before forcing a decay, the generator checks `m_parent > 2*m_chi`.  Closed combinations write a zero summary row with `stop_reason = kinematically_closed` and do not poison Pythia with impossible decay tables.

## Build

```bash
cmake -S . -B build
cmake --build build -j
```

or directly:

```bash
g++ -O2 -std=c++17 -pthread allInOne_multimeson.cc -o allInOne_multimeson \
  $(root-config --cflags --libs) \
  $(pythia8-config --cxxflags --libs)
```

## Run one local job

Fixed generated events are preferred for normalization:

```bash
./allInOne_multimeson 1000 4 10000 0.01 2x2 pi0 light_mesons outputs/pi0_m0p01.root \
  --mode fixed-events --n-events 10000 --write-spectra accepted
```

Target-accepted mode remains available for convenience:

```bash
./allInOne_multimeson 1000 4 100 0.45 2x2 phi light_mesons outputs/phi_m0p45.root \
  --mode target-accepted --n-accepted 100 --max-events 1000000
```

Useful options:

- `--write-spectra all|accepted|none`
- `--spectra-prescale N`
- `--beam-config beam.config`
- `--momentum-config momentum.config`
- `--production-config configs/charmonium.config`
- `--print-pythia-settings`
- `--quiet` / `--batch`
### Production config diagnostics

`production_mode=charmonium` explicitly disables SoftQCD/minimum-bias switches after reading the common beam setup and before enabling the charmonium card.  The default editable card is:

```text
configs/charmonium.config
```

The default charmonium settings used by that card are:

```text
HardQCD:hardccbar = on
Charmonium:all = on
```

You can test alternate charmonium steering without recompiling:

```bash
./allInOne_multimeson 1000 1 10000 1.00 2x2 jpsi charmonium outputs/test_jpsi_1gev.root \
  --mode fixed-events --n-events 10000 --write-spectra accepted \
  --production-config configs/charmonium.config --print-pythia-settings
```

`--print-pythia-settings` prints Pythia's changed settings before `init()`, which is useful for verifying that `SoftQCD:all = off` in charmonium mode.

The `n_pythia_next_failures` counter is incremented whenever `pythia.next()` returns `false`.  Pythia warning messages that do not cause a false return are not included in this counter.


Intentional breaking change: the old `parentList=all` multi-emitter CLI is replaced by a single `emitterName productionMode` per job.  The executable name `allInOne_multimeson` is retained for backwards build compatibility; CMake also provides an `allInOne_emitter` alias target.

## ROOT output

Each job writes one ROOT file with two TTrees.

### `mcp_summary`

Important branches include:

- metadata: `run_id`, `job_id`, `thread_id`, `seed`, `mcp_mass_GeV`, `emitter_pdg`, `emitter_name`, `emitter_type`, `production_mode`, `production_mode_name`, `geometry_id`, `geometry_name`, `is_kinematically_open`
- counters: `n_events_generated`, `n_pythia_next_failures`, `n_emitter_total`, `n_mcp_total`, `n_mcp_accepted`, `n_mcp_wrong_mother`, `n_mcp_no_mother`
- derived: `emitter_per_event`, `acceptance_fraction`, `acceptance_uncertainty_binomial`, `stop_reason`, `stopping_mode`, `spectra_prescale`

Backward-compatible aliases such as `mcp_mass`, `parent_pdg`, `parent_type`, `parent_name`, `n_parent_total`, and `parent_yield_per_event` are also written.

### `mcp_spectra`

One row per retained MCP according to `--write-spectra` and `--spectra-prescale`, including event metadata, MCP/mother identity, lab four-vectors, angular variables, detector projection, mother four-vector, and feed-down diagnostics.

Quick ROOT diagnostic command:

```bash
root -l -b -q -e '
TFile f("outputs/test_jpsi_1gev.root");
auto s = (TTree*)f.Get("mcp_summary");
auto k = (TTree*)f.Get("mcp_spectra");
s->Scan("mcp_mass_GeV:emitter_name:production_mode_name:is_kinematically_open:n_events_generated:n_emitter_total:n_mcp_total:n_mcp_accepted:n_mcp_wrong_mother:n_mcp_no_mother:n_pythia_next_failures:stop_reason:stopping_mode_name", "", "colsize=18", 20);
std::cout << "spectra entries = " << (k ? k->GetEntries() : -1) << std::endl;
'
```


## NERSC-style workflow

Generate a job table:

```bash
scripts/make_job_list.py --masses masses.txt --shards 10 --n-events 100000 --geometry 2x2 --output job_list.csv
```

Run one row locally or from SLURM:

```bash
JOB_LIST=job_list.csv EXEC=./allInOne_multimeson NTHREADS=1 scripts/run_one_job.sh 0
```

Submit an array on Perlmutter CPU after editing account/array size placeholders:

```bash
mkdir -p logs/slurm
sbatch --array=0-999 scripts/submit_nersc_array.slurm
```

## Aggregation

```bash
scripts/aggregate_outputs.py outputs/raw --csv outputs/aggregate_summary.csv
```

The aggregator groups by `mcp_mass_GeV`, `emitter_pdg`, `production_mode`, and `geometry_id`, keeping light-meson and charmonium normalization separate.

## Plotting

The plotting macro consumes raw ROOT files containing `mcp_summary` (or ROOT aggregates with the same tree name):

```cpp
.x plot_acceptance_vs_mass_multimeson.C
plot_acceptance_vs_mass_multimeson("outputs/raw/emitter_pi0/mass_0p010000/job_000000.root", 1.5e19, 1e-2);
```

It plots acceptance vs mass and an accepted flux proxy.  Placeholder weights are applied only at plotting time:

- pseudoscalars: `(1 - 4*mchi^2/mM^2)^3`
- vectors: `sqrt(1 - 4*mchi^2/mV^2)`

Absolute normalizations are intentionally editable placeholders, not PDG-perfect results.

## Validation checklist

Suggested smoke tests:

1. `pi0`, `mchi=0.01 GeV`: expect MCPs.
2. `pi0`, `mchi=0.10 GeV`: expect closed and skipped.
3. `phi`, `mchi=0.45 GeV`: expect open, likely low statistics.
4. `phi`, `mchi=0.60 GeV`: expect closed.
5. `jpsi`, `mchi=1.0 GeV`: expect open if local charmonium production is configured correctly.  Smoke-test command:

   ```bash
   ./allInOne_multimeson 1000 1 10000 1.00 2x2 jpsi charmonium outputs/test_jpsi_1gev.root \
     --mode fixed-events --n-events 10000 --write-spectra accepted --quiet
   ```

   Expected: `is_kinematically_open = 1`, `n_events_generated = 10000`, clearly nonzero `n_emitter_total`, `n_mcp_total = 2 * n_emitter_total` unless Pythia event-record status semantics require further documented filtering, and no `should not combine softQCD processes with hard ones` warning.
6. In a `pi0` run, `n_mcp_wrong_mother` should be zero; MCP mothers should be PDG 111, including feed-down pi0s.

## Caveats

- MCP branching ratios are forced to 1 in generation for statistics.
- Epsilon scaling, physical branching ratios, and detector response are applied later.
- Charmonium production is isolated for bookkeeping; validate its Pythia tune before absolute interpretation.
- DY is reserved but not generated in this revision.
