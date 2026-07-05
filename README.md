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
- `--diagnostics`
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



## Why charmonium is normalized separately

`light_mesons` mode and `charmonium` mode are intentionally **different production samples**.

`light_mesons` mode is a SoftQCD/minimum-bias-style sample. In this repository's current convention, it is used to estimate light neutral meson production per generic proton interaction. This is the same convention used in the original pi0-only pipeline, where a meson yield per generated event was scaled by POT.

`charmonium` mode is not a generic proton-interaction sample. It explicitly enables hard charm/charmonium production, for example the default card:

```text
HardQCD:hardccbar = on
Charmonium:all = on
```

This biases generation toward charm/charmonium events so that the J/psi acceptance and spectra can be measured with useful statistics. That is good for efficiency, but it means raw event counts in `charmonium` mode cannot be added directly to raw event counts in `light_mesons` mode.

Consequences:

- `n_events_generated` in `charmonium` mode is not directly comparable to `n_events_generated` in `light_mesons` mode.
- Do not add light-meson and charmonium yields by raw generated-event counts.
- Charmonium should be normalized using a generated/known production cross section or an external production model.
- The generator records `production_mode`, `sigma_gen_mb`, `sigma_err_mb`, and `weight_per_event_mb` to support later normalization.
- The absolute J/psi normalization is therefore **Pythia-card/tune dependent** and should be validated before using it as a publication-quality prediction.

The normalized plotting layer treats charmonium as an importance-sampled component. The scale factor used to convert a charmonium-biased event sample into the light-meson/minimum-bias convention is

$$
S_{J/\psi}
=
\frac{\langle \sigma_{\rm charmonium} \rangle}{\langle \sigma_{\rm SoftQCD} \rangle}.
$$

Here the angle brackets mean the **mean generated cross section per input file/shard**, not the sum over files. This distinction is important: `sigma_gen_mb_sum` grows if you run more shards, while the physical process cross section should not.

For an aggregate containing equivalent shards,

$$
\langle \sigma_{\rm mode} \rangle
=
\frac{\sum_i \sigma_{{\rm gen},i}}{N_{\rm input\ files}}.
$$

Then the default process scales are

$$
S_{\rm light\ mesons}=1,
\qquad
S_{\rm charmonium}=\frac{\langle \sigma_{\rm charmonium}\rangle}{\langle \sigma_{\rm SoftQCD}\rangle}.
$$

A useful sanity check is that doubling the number of equivalent shards should roughly double `n_events_generated`, `n_mcp_total`, and `n_mcp_accepted`, but it should **not** double the mean cross section used for normalization.

Intentional breaking change: the old `parentList=all` multi-emitter CLI is replaced by a single `emitterName productionMode` per job.  The executable name `allInOne_multimeson` is retained for backwards build compatibility; CMake also provides an `allInOne_emitter` alias target.

## ROOT output

Each job writes one ROOT file with two TTrees.

### `mcp_summary`

Important branches include:

- metadata: `run_id`, `job_id`, `thread_id`, `seed`, `mcp_mass_GeV`, `emitter_pdg`, `emitter_name`, `emitter_type`, `production_mode`, `production_mode_name`, `geometry_id`, `geometry_name`, `is_kinematically_open`
- counters: `n_events_generated`, `n_pythia_next_failures`, `n_emitter_record_entries`, `n_emitter_decayed_to_mcp`, `n_mcp_pairs`, `n_mcp_pairing_anomalies`, `n_emitter_total`, `n_mcp_total`, `n_mcp_accepted`, `n_mcp_wrong_mother`, `n_mcp_no_mother`
- derived/cross-section metadata: `emitter_per_event`, `acceptance_fraction`, `acceptance_uncertainty_binomial`, `sigma_gen_mb`, `sigma_err_mb`, `weight_per_event_mb`, `stop_reason`, `stopping_mode`, `spectra_prescale`

Backward-compatible aliases such as `mcp_mass`, `parent_pdg`, `parent_type`, `parent_name`, `n_parent_total`, and `parent_yield_per_event` are also written.  `n_emitter_total` / `n_parent_total` now alias `n_emitter_decayed_to_mcp`, which is the more physical forced-decay count.



### Event-record emitter bookkeeping

Pythia event records can contain several entries with the same PDG ID for one physical object because of history/status copies.  Therefore:

- `n_emitter_record_entries` counts all event-record entries with `abs(id) == emitter_pdg`; it is useful for status diagnostics but is not necessarily a physical decay count.
- `n_emitter_decayed_to_mcp` counts unique MCP decay mothers using the event-local mother together with the thread/event identity, effectively `(thread_id, event_index, mother_index)`.
- `n_mcp_pairs` counts unique MCP mothers with exactly two MCP daughters.
- `n_mcp_pairing_anomalies` counts unique MCP mothers with a daughter count different from two.
- For these forced two-body vector decays, `n_mcp_total` should usually equal `2 * n_emitter_decayed_to_mcp`; deviations print a warning and are stored in the anomaly counter.

The file also contains an `emitter_status_counts` tree with `mcp_mass_GeV`, `emitter_pdg`, `production_mode`, `geometry_id`, `status_code`, and `count` to diagnose Pythia status/history copies.

### `mcp_spectra`

One row per retained MCP according to `--write-spectra` and `--spectra-prescale`, including event metadata, MCP/mother identity, lab four-vectors, angular variables, detector projection, mother four-vector, and feed-down diagnostics.

Quick ROOT diagnostic command:

```bash
root -l -b -q -e '
TFile f("outputs/test_jpsi_1gev.root");
auto s = (TTree*)f.Get("mcp_summary");
auto k = (TTree*)f.Get("mcp_spectra");
s->Scan("mcp_mass_GeV:emitter_name:production_mode_name:is_kinematically_open:n_events_generated:n_emitter_record_entries:n_emitter_decayed_to_mcp:n_mcp_pairs:n_mcp_pairing_anomalies:n_mcp_total:n_mcp_accepted:n_pythia_next_failures:sigma_gen_mb:stop_reason:stopping_mode_name", "", "colsize=18", 20);
std::cout << "spectra entries = " << (k ? k->GetEntries() : -1) << std::endl;
'
```

Unique MCP-mother diagnostic using `(event_index,mother_index)` for a single-thread job:

```bash
root -l -b -q -e '
#include <set>
#include <utility>
TFile f("outputs/test_jpsi_1gev.root");
auto k = (TTree*)f.Get("mcp_spectra");
std::set<std::pair<int,int>> mothers;
int event_index = 0, mother_index = 0, mother_pdg = 0;
k->SetBranchAddress("event_index", &event_index);
k->SetBranchAddress("mother_index", &mother_index);
k->SetBranchAddress("mother_pdg", &mother_pdg);
for (Long64_t i = 0; i < k->GetEntries(); ++i) { k->GetEntry(i); if (mother_pdg == 443) mothers.insert({event_index, mother_index}); }
std::cout << "unique J/psi MCP mothers = " << mothers.size() << std::endl;
'
```

Run with `--diagnostics` to print retained-spectra geometry summaries, including x/y detector projection bands and min/max/mean/RMS for `theta_x_rad`, `theta_y_rad`, `x_at_detector_m`, and `y_at_detector_m`.

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

The aggregator groups by `mcp_mass_GeV`, `emitter_pdg`, `production_mode`, and `geometry_id`, keeping light-meson and charmonium normalization separate.  It preserves the record-entry/decayed-emitter/pairing counters and carries cross-section metadata to CSV/ROOT for diagnostics; do not interpret summed charmonium cross sections as an absolute normalization without validation.

## Plotting

The plotting macros consume either raw ROOT files containing `mcp_summary` or aggregate ROOT files with the same tree name.

Legacy/proxy plotting:

```cpp
.x plot_acceptance_vs_mass_multimeson.C
plot_acceptance_vs_mass_multimeson("outputs/raw/emitter_pi0/mass_0p010000/job_000000.root", 1.5e19, 1e-2);
```

The legacy plot is useful for quick debugging, but it should be interpreted as a **proxy** unless the normalization choices below are applied.

Normalized plotting:

```cpp
.x plot_normalized_mcp_yield.C
plot_normalized_mcp_yield("outputs/aggregate_summary_4h.root", 1.5e19, 1e-2);
```

This produces:

1. acceptance vs MCP mass per emitter;
2. accepted MCP yield vs MCP mass per emitter;
3. total accepted MCP yield vs MCP mass, summed over emitters;
4. labels showing the chosen `N_POT`, epsilon, and the fact that the branching constants are configurable model inputs.

The normalized-yield plot should be read as a pipeline-level physics prediction using the current reference constants and phase-space approximations. It is not yet a final sensitivity result until the branching constants, production model, detector response, and background treatment are validated.

## Validation checklist

Suggested smoke tests:

1. `pi0`, `mchi=0.01 GeV`: expect MCPs.
2. `pi0`, `mchi=0.10 GeV`: expect closed and skipped.
3. `phi`, `mchi=0.45 GeV`: expect open, likely low statistics.
4. `phi`, `mchi=0.60 GeV`: expect closed.
5. `jpsi`, `mchi=1.0 GeV`: expect open if local charmonium production is configured correctly.  Smoke-test command:

   ```bash
   ./allInOne_multimeson 1000 1 10000 1.00 2x2 jpsi charmonium outputs/test_jpsi_1gev.root \
     --mode fixed-events --n-events 10000 --write-spectra all --spectra-prescale 1 --diagnostics --quiet
   ```

   Expected: `is_kinematically_open = 1`, `n_events_generated = 10000`, clearly nonzero `n_emitter_decayed_to_mcp`, `n_mcp_total = 2 * n_emitter_decayed_to_mcp`, `n_mcp_pairing_anomalies = 0`, and no `should not combine softQCD processes with hard ones` warning. `n_emitter_record_entries` may be larger than `n_emitter_decayed_to_mcp` because of Pythia event-record/status copies.
6. In a `pi0` run, `n_mcp_wrong_mother` should be zero; MCP mothers should be PDG 111, including feed-down pi0s.

## Caveats

- MCP branching ratios are forced to 1 in generation for statistics.
- Epsilon scaling, physical branching ratios, and detector response are applied later.
- Charmonium production is isolated for bookkeeping; validate its Pythia tune before absolute interpretation.
- DY is reserved but not generated in this revision.

## Normalization to POT

The normalized plotting and toy-MC export tools use explicit exotic branching weights. The default numbers are **model inputs**, not fit parameters. They are intended to make the pipeline reproducible and easy to update as the physics model is refined.

There are three different kinds of normalization ingredients:

1. **Reference branching ratios**: fixed lookup values, normally taken from PDG or another documented source. These are not intended to float in the analysis, but they should be easy to update.
2. **Phase-space factors**: analytic/model approximations for the MCP-mass dependence. These are not arbitrary, but the current formulas are approximations that may later be replaced by a more exact integral.
3. **Production cross sections**: generator/model dependent. This matters especially for charmonium, which is generated with a biased hard/charmonium process card.

For a parent/emitter and MCP mass, the accepted MCP yield is computed as

$$
N_{\rm acc}(P,m_\chi,\varepsilon)
=
N_{\rm POT}
\times
S_P
\times
\frac{N_{\chi,{\rm accepted}}}{N_{\rm events}}
\times
{\rm Br}_{\rm exotic}(P,m_\chi,\varepsilon).
$$

Here:

- $P$ is the selected emitter species;
- $S_P$ is the production-mode scale factor;
- $N_{\chi,{\rm accepted}}$ is `n_mcp_accepted`;
- $N_{\rm events}$ is `n_events_generated`.

`n_mcp_accepted` already counts MCP particles. Therefore this formula does **not** multiply by two again.

For `light_mesons` mode, the SoftQCD/minimum-bias-style sample is treated as generic proton interactions under the same convention used in the earlier acceptance plots:

$$
S_{\rm light\ mesons}=1.
$$

For `charmonium` mode, the sample is a biased hard/charmonium sample. Its raw event count cannot be added directly to light-meson event counts. The plotting/export tools use

$$
S_{\rm charmonium}
=
\frac{\langle \sigma_{\rm charmonium} \rangle}
       {\langle \sigma_{\rm SoftQCD} \rangle}.
$$

The mean cross section is computed per input file/shard:

$$
\langle \sigma_{\rm mode} \rangle
=
\frac{\sigma_{\rm gen,mb,sum}}{N_{\rm input\ files}}.
$$

Do **not** use `sigma_gen_mb_sum` directly as a process cross section. It grows with the number of shards. The aggregator stores both sum and mean values so the normalization remains stable when more shards are added.

### Exotic branching weights

Let

$$
r = \frac{m_\chi^2}{m_P^2}.
$$

The current default pseudoscalar approximation for

$$
P \to \gamma\chi\bar{\chi}
$$

is

$$
{\rm Br}_{\rm exotic}^{\rm pseudo}(P,m_\chi,\varepsilon)
=
\varepsilon^2\,\alpha\,{\rm Br}_{\rm ref}(P)
\left(1-4r\right)^3,
\qquad
2m_\chi < m_P.
$$

For closed channels, $2m_\chi \geq m_P$, the branching weight is set to zero.

The current default vector approximation for

$$
V \to \gamma^* \to \chi\bar{\chi}
$$

is

$$
{\rm Br}_{\rm exotic}^{\rm vector}(V,m_\chi,\varepsilon)
=
\varepsilon^2\,\alpha\,{\rm Br}_{\rm ref}(V)
\,\beta\,(1+2r),
\qquad
\beta = \sqrt{1-4r},
\qquad
2m_\chi < m_V.
$$

Again, closed channels receive zero weight.

The older pseudoscalar integral approximation can be kept as an optional future switch, but the default power-law form above preserves continuity with the existing pipeline.

### Default reference constants

The defaults below are the constants used by the current plotting code. They should be documented in code comments and validated against the desired PDG release before final results.

| Emitter | PDG | Type | Mass [GeV] | Reference branching normalization |
|---|---:|---|---:|---:|
| $\pi^0$ | 111 | pseudoscalar | 0.1349770 | ${\rm Br}(\pi^0\to\gamma\gamma)=0.98823$ |
| $\eta$ | 221 | pseudoscalar | 0.5478620 | ${\rm Br}(\eta\to\gamma\gamma)=0.39410$ |
| $\eta'$ | 331 | pseudoscalar | 0.9577800 | ${\rm Br}(\eta'\to\gamma\gamma)=0.02220$ |
| $\rho^0$ | 113 | vector | 0.7752600 | ${\rm Br}(\rho^0\to e^+e^-)=4.72\times 10^{-5}$ |
| $\omega$ | 223 | vector | 0.7826500 | ${\rm Br}(\omega\to e^+e^-)=7.36\times 10^{-5}$ |
| $\phi$ | 333 | vector | 1.0194610 | ${\rm Br}(\phi\to e^+e^-)=2.973\times 10^{-4}$ |
| $J/\psi$ | 443 | vector | 3.0969000 | ${\rm Br}(J/\psi\to e^+e^-)=5.971\times 10^{-2}$ |

The pseudoscalar constants use the two-photon branching fraction as the reference electromagnetic normalization. The vector constants use the dielectron branching fraction as the reference electromagnetic normalization. This is a compact phenomenological model choice; the code should keep these constants centralized so a more exact model can replace them later.

Run the normalized-yield macro with:

```cpp
.x plot_normalized_mcp_yield.C
plot_normalized_mcp_yield("outputs/aggregate_summary_4h.root", 1.5e19, 1e-2);
```

It produces acceptance vs mass, accepted MCP yield per emitter, and total accepted MCP yield summed over emitters, with labels noting $N_{\rm POT}$, $\varepsilon$, and the placeholder branching normalization.

## Spectra and angular plots

Use the spectra plotting helper on raw ROOT files containing `mcp_spectra`:

```bash
scripts/plot_spectra.py \
  --spectra-dir outputs/raw_4h \
  --output-dir plots/spectra_4h \
  --accepted-only
```

The script makes `theta_x_rad` vs `theta_y_rad`, theta, `pz_GeV` vs theta, energy, detector-plane x/y, and accepted-only x/y hit-map plots.  The detector-plane plot draws the 2x2 boxes:

```text
left:  x in [-0.65, -0.05], y in [-0.7, 0.7]
right: x in [ 0.05,  0.65], y in [-0.7, 0.7]
```

Production runs with `WRITE_SPECTRA=accepted` are immediately useful for hit maps.  Full angular distributions require spectra-focused runs such as `--write-spectra all --spectra-prescale N`.

## Toy MC spectra export

Export weighted MCP rows for toy detector MC with:

```bash
scripts/export_toymc_spectra.py \
  --summary outputs/aggregate_summary_4h.csv \
  --spectra-dir outputs/raw_4h \
  --pot 1.5e19 \
  --epsilon 1e-2 \
  --output-csv outputs/toymc_mcp_spectra_eps1e-2.csv \
  --accepted-only
```

The output CSV includes kinematics, detector projection, `process_scale`, `br_exotic`, `spectra_prescale`, and weights. For each retained MCP spectrum row,

$$
w_{\rm event}
=
N_{\rm POT}
\times
S_P
\times
{\rm Br}_{\rm exotic}(P,m_\chi,\varepsilon)
\times
\frac{f_{\rm prescale}}{N_{\rm events}(P,m_\chi,{\rm mode})}.
$$

The per-POT weight is

$$
w_{\rm per\ POT}=\frac{w_{\rm event}}{N_{\rm POT}}.
$$

For accepted-only spectra with `spectra_prescale = 1`, summing `event_weight` over rows with `passed_geometry = 1` should reproduce the normalized accepted-yield plot for that mass/emitter/mode. For all-produced spectra with `spectra_prescale > 1`, the prescale factor accounts for downsampling.

The export script validates this by grouping rows with

$$
(m_\chi,\;{\rm emitter\_pdg},\;{\rm production\_mode},\;{\rm geometry\_id}).
$$

The production mode is always part of the key so light-meson and charmonium normalizations are not mixed.

## Current caveats and interpretation notes

- MCP branching ratios are forced to 1 in generation for statistics. Physical branching weights are applied only in post-processing.
- Epsilon scaling, physical branching ratios, detector response, and background treatment are not handled by the generator.
- The default branching constants are documented reference inputs, not fit parameters, but they still need validation against the desired PDG release/model choice.
- The phase-space functions are analytic approximations. They are useful for a consistent pipeline but may be replaced by a more exact model.
- The J/psi/charmonium component is generated with a biased hard/charmonium sample and must be normalized using cross sections or an external production model. Do not compare raw `n_events_generated` between `charmonium` and `light_mesons`.
- `sigma_gen_mb_sum` is diagnostic. Use `sigma_gen_mb_mean` or an equivalent per-file mean for process normalization.
- DY is reserved but not generated in this revision.
- A final sensitivity curve still needs detector response, exposure treatment, backgrounds, and systematic uncertainties.
