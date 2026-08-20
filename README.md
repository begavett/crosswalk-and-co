# crosswalk-and-co

A simulation study comparing **scale cocalibration** and **score crosswalking** as methods for harmonizing cognitive test scores across groups — developed at PsyMCA 2025.

**Authors:** Brandon Gavett, Emma Nichols, Nathan Ryder
**Collaborators:** Hoda Abdel Magid, Sarah Forrester, Xiaqing Jiang, Katrina Kezios, Adina Zeki Al Hazzouri, Alden Gross, Jason Gantenberg, Sarah Ackley

---

## Overview

When cognitive assessments are administered to different groups (e.g., studies using different test batteries, or populations differing by language or culture), the resulting scores may not be directly comparable. Two common approaches to address this are:

- **Cocalibration** — placing scores from different groups onto a common IRT scale using linking (anchor) items.
- **Effect size crosswalking** — statistically mapping effect sizes from one test to another using the [`cogxwalkr`](https://github.com/jrgant/cogxwalkr) package.

This project simulates cognitive item response data under four scenarios varying the degree of differential item functioning (DIF) and anchor item strength, then benchmarks each method's ability to recover a known education–cognition regression coefficient.

### Simulation scenarios

| Scenario | DIF    | Anchor strength |
|----------|--------|-----------------|
| 1        | None   | Strong          |
| 2        | None   | Weak            |
| 3        | Weak   | Strong          |
| 4        | Strong | Strong          |

### Simulated groups

Each replication simulates three groups of size `N` from the same latent education→cognition model:

| Group   | Item parameters                                   | Role |
|---------|---------------------------------------------------|------|
| Group 1 | `b_m` sheet of `data/b_m.xlsx`                     | Reference / focal group for cocalibration |
| Group 2 | `b_m_group2` sheet (DIF-perturbed parameters)      | Focal / reference group for cocalibration |
| Group 3 | Item-by-item random mixture of Group 1 and Group 2 | **Crosswalk (linking) sample** — the external sample in which the crosswalk between measures is estimated |

Cocalibration uses Groups 1 and 2. The `cogxwalkr` arm estimates its crosswalk in Group 3 and applies it to the education–cognition association estimated in Groups 1 and 2.

### Estimands compared

Each replication returns, per scenario and target group:

- `coef` — the harmonized education coefficient (with `lwr`/`upr` when bootstrap CIs are enabled)
- `truecoefs` — the target coefficient (education→score regression on the target group's own IRT score)
- `naivecoefs` — the unharmonized coefficient (using whatever outcome was observed, no harmonization)
- `b1` — the true simulated effect size

Bias is therefore evaluated two ways: against the **simulated** truth (`coef - b1`) and against the **target** coefficient (`coef - truecoefs`).

---

## Repository structure

```
crosswalk-and-co/
├── sim_cwxco.R              # Main simulation driver (parallel, local/laptop)
├── simstudy_example.sh      # SLURM array-job template for cluster runs
├── getresults_example.R     # Collects per-replication RDS files into one results object
├── summarize-results.R      # Exploratory bias boxplots/barplots for a single results file
├── paperplots.R             # Publication figures (bias, mean bias, CI coverage)
├── crosswalk.R              # Standalone/earlier crosswalk-only script (superseded by src/helper.R)
├── src/
│   ├── helper.R             # Core functions: simCog, cocalibrate, getscore, getcrosswalk, sim_cwxco, ...
│   └── crosswalk_helper.R   # Earlier standalone getscore/getcrosswalk wrappers
├── data/
│   ├── b_m.xlsx             # IRT item parameters (sheets: b_m, b_m_group2)
│   ├── scenarios.csv        # Scenario × group → item map (4 scenarios × 2 groups)
│   ├── RecodedItemName_m.csv# Item name lookup (40 names; 31 used after dropping `rav*b` items)
│   └── a_m.csv              # Discrimination parameters (reference copy)
├── models/
│   ├── dem_logr.Rds         # Logistic model for dementia probability (includes education)
│   ├── dem_logr_noedu.Rds   # Same model without education
│   ├── cc_models.Rds        # Saved cocalibration model objects (used for test information plots)
│   ├── allresn*_d*_*.RDS    # Aggregated simulation results (see naming below)
│   └── archive/             # Earlier result sets, kept for reference
└── plots/                   # Generated figures (jpg/pdf)
```

### Results file naming

Aggregated result files in `models/` follow `allresn{N}_d{b1×10}_{variant}[_{month}].RDS`:

- `n500`, `n1000`, `n5000` — sample size per group
- `d2`, `d4` — true effect size `b1` of 0.2 or 0.4
- `ci` — run with bootstrap CIs for the `cogxwalkr` method; `noci` — without
- `sum` — variant using sum scores as the crosswalk input
- `jun2026` / `aug2026` — run date; `withAug25package` marks a run pinned to an older `cogxwalkr` build

`paperplots.R` currently reads the six files matching `allresn[0-9]*_d[0-9]_ci_aug2026`.

---

## Dependencies

Requires **R ≥ 4.3**. Install packages with:

```r
install.packages("pacman")
pacman::p_load(
  dplyr, tidyr, magrittr, purrr, psych, mirt, data.table, ggplot2,
  readxl, readr, stringr, simDAG, cogxwalkr, forcats,
  future, future.apply, progressr, patchwork, ggmirt, parameters,
  ggtext, binom
)
```

`cogxwalkr` is installed from GitHub:

```r
remotes::install_github("jrgant/cogxwalkr@04a4cb7725736ca7ad9c860f69c9277f0095b31a")
```

`ggmirt` (used only for the test information plots at the end of `sim_cwxco.R`) is also from GitHub:

```r
remotes::install_github("masurp/ggmirt")
```

---

## Usage

### 1. Set your working directory

`sim_cwxco.R`, `summarize-results.R`, and `paperplots.R` each detect the current user and set `code_dir` accordingly. Add your own entry:

```r
if (user == "yourname") {
  code_dir <- "/path/to/crosswalk-and-co/"
}
```

### 2. Run locally (parallel)

`sim_cwxco.R` runs replications with `future_Map()`. Defaults are 1,000 iterations at N = 1,002 per group across 8 workers:

```r
iter_start   <- 1
niter        <- 1000
sample_size  <- 1002
plan(multisession, workers = 8)
```

A single serial run for testing:

```r
sim_cwxco(iter = 1, N = 1001, prop_high_edu = .30, b0 = 0, b1 = .2, dem_prob_cut = .9)
```

Results are collected with:

```r
results_df <- rbindlist(results, idcol = "rep") %>% mutate(n_sample = sample_size)
```

If `models/cc_models.Rds` is present, the tail of the script draws test information curves for all four scenarios via `ggmirt::testInfoPlot()`.

### 3. Run on a cluster (SLURM)

`simstudy_example.sh` is a template array job that runs one replication per array task and writes `result{task}.RDS` into a per-condition output directory:

```bash
sbatch simstudy_example.sh
```

Edit `n_sample`, `b1_val`, `CI`, `seedpad`, and the `code_dir` / `outdir` paths at the top of the heredoc. The array index becomes the seed, so `seedpad` lets you extend a run without reusing seeds.

Then combine the per-task files with `getresults_example.R`, which reads all `result{rep}.RDS` files for a given N and effect size, reports any missing replications, and saves a single aggregated RDS.

### 4. Summarize and plot

- `summarize-results.R` — `summarize_sims(outfilename, ref_grp, plot_type = c("boxplot", "barplot"), vs = c("truth", "target"))` for quick looks at one results file.
- `paperplots.R` — builds the manuscript figures: bias distributions, percent-bias means with t-based CIs, an education-condition comparison, and bootstrap CI coverage. Figures are written to `plots/` with the run date in the filename.

---

## Key functions (`src/helper.R`)

| Function | Description |
|----------|-------------|
| `simCog()` | Simulates IRT-based cognitive item responses with a latent regression on education, then returns true thetas plus estimated factor scores. |
| `cocalibrate()` | Places two groups on a common IRT scale by estimating the reference group first, then fixing linking-item parameters in the focal group while freeing its mean and variance. |
| `getscore()` | Fits a unidimensional graded response model and appends factor scores. |
| `getsumscore()` | Appends a simple sum score across a set of items (used as the crosswalk input). |
| `getcrosswalk()` | Wraps `cogxwalkr::crosswalk()`, conditioning on dementia status; optionally bootstraps (500 resamples) when `CI = TRUE`. |
| `cocalibrated_regressions()` | Fits education–cognition regressions on cocalibrated scores and extracts coefficients and CIs. |
| `get_scenario_data()` | Extracts scenario-specific item subsets plus key variables from the full dataset. |
| `sim_cwxco()` | Master function: simulates all three groups, runs cocalibration and crosswalks across the four scenarios, and returns combined results. |

---

## Simulation parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `iter` | 1 | Replication index; also used as the random seed (`set.seed(iter)`) and bootstrap seed |
| `N` | 5002 | Sample size per group |
| `prop_high_edu` | 0.30 | Proportion with high education |
| `b0` | 0 | Intercept for the latent education regression |
| `b1` | 0.20 | Slope (true education effect on cognition) |
| `dem_prob_cut` | 0.90 | Probability threshold for assigning dementia status |
| `CI` | `FALSE` | Bootstrap CIs (500 resamples) for the `cogxwalkr` method; substantially increases runtime |

Dementia status is assigned two ways in each run — from `dem_logr` (education included) and `dem_logr_noedu` — so the crosswalk can be conditioned either way. The resulting rows are flagged by the `demwithedu` column.

---

## Notes

- Parallel execution uses the [`future`](https://future.futureverse.org/) framework. On HPC, prefer the SLURM array approach in `simstudy_example.sh` over `multisession`.
- `future.globals.maxSize` is set to 1,500 MB in `sim_cwxco.R`; raise it if you hit globals-size errors.
- Memory and runtime scale with `N × niter`, and `CI = TRUE` multiplies crosswalk cost by the number of bootstrap resamples. Reduce both for testing.
- Seeds are set per iteration for reproducibility.
- `plots/` is listed in `.gitignore`, so newly generated figures are not tracked by default; the existing files in that directory were committed before the rule was added.

## Known rough edges

These are inconsistencies in the current code, not features:

- `sim_cwxco.R` passes `dem_on_edu = FALSE` to `sim_cwxco()`, but that argument is no longer in the function signature (both dementia definitions are now computed unconditionally). The call will error until the argument is dropped.
- `simstudy_example.sh` sources `src/helper_sum.R`, which is not in this repository. Point it at `src/helper.R`, or add the sum-score helper file used to produce the `*_sum.RDS` results.
- `summarize-results.R` looks for results in `data/`, but aggregated results now live in `models/`. It also ends with a top-level call to `summarize_sims("allresn5000_d1.Rds", ...)`, a filename that no longer exists — the script will fail if sourced as-is.
- The `cc_models.Rds` save block at the end of `sim_cwxco()` is commented out, so the plotting section of `sim_cwxco.R` only runs if the file already exists.
- `crosswalk.R` and `src/crosswalk_helper.R` define earlier versions of `getscore()`/`getcrosswalk()` that are superseded by `src/helper.R` and depend on objects (`item_lists`, `combos`, `scenario_data`) that only exist inside `sim_cwxco()`. They are kept for reference and are not part of the current pipeline.
