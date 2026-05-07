# crosswalk-and-co

A simulation study comparing **scale cocalibration** and **score crosswalking** as methods for harmonizing cognitive test scores across groups — developed at PsyMCA 2025.

**Authors:** Brandon Gavett, Emma Nichols, Nathan Ryder  
**Collaborators:** Hoda Abdel Magid, Sarah Forrester, Xiaqing Jiang, Katrina Kezios, Adina Zeki Al Hazzouri

---

## Overview

When cognitive assessments are administered to different groups (e.g., studies using different test batteries, or populations differing by language or culture), the resulting scores may not be directly comparable. Two common approaches to address this are:

- **Cocalibration** — placing scores from different groups onto a common IRT scale using linking (anchor) items.
- **Score crosswalking** — statistically mapping scores from one test to another using the [`cogxwalkr`](https://github.com/jrgant/cogxwalkr) package.

This project simulates cognitive item response data under four scenarios varying the degree of differential item functioning (DIF) and anchor item strength, then benchmarks each method's ability to recover a known education–cognition regression coefficient.

### Simulation Scenarios

| Scenario | DIF      | Anchor Strength |
|----------|----------|-----------------|
| 1        | None     | Strong          |
| 2        | None     | Weak            |
| 3        | Weak     | Strong          |
| 4        | Strong   | Strong          |

---

## Repository Structure

```
crosswalk-and-co/
├── crosswalk.R          # Standalone crosswalk analysis using cogxwalkr
├── sim_cwxco.R          # Main simulation entry point (parallel execution)
├── src/
│   ├── helper.R         # Core functions: simCog, cocalibrate, sim_cwxco, etc.
│   └── crosswalk_helper.R  # Utility wrappers for getscore and getcrosswalk
├── data/
│   ├── b_m.xlsx         # IRT item parameters for Groups 1 & 2
│   ├── scenarios.csv    # Scenario-to-item mapping
│   └── RecodedItemName_m.csv  # Item name lookup table
├── models/
│   ├── dem_logr.Rds     # Logistic regression model for dementia probability
│   ├── dem_logr_noedu.Rds
│   └── cc_models.Rds    # Saved cocalibration model objects (generated on first run)
└── summarize-results.R  # (Optional) Result summarization script
```

---

## Dependencies

Requires **R ≥ 4.3**. Install packages with:

```r
install.packages("pacman")
pacman::p_load(
  dplyr, tidyr, magrittr, psych, mirt, data.table, ggplot2,
  readxl, readr, stringr, simDAG, cogxwalkr, forcats,
  future, future.apply, progressr, patchwork, ggmirt, parameters
)
```

The `cogxwalkr` package is installed from GitHub:

```r
remotes::install_github("jrgant/cogxwalkr")
```

---

## Usage

### 1. Configure your working directory

At the top of `sim_cwxco.R`, set your `code_dir` path. The script detects the current user automatically for known users (`brandon`, `emmanich`); add your own entry if needed:

```r
if (user == "yourname") {
  code_dir <- "/path/to/crosswalk-and-co/"
}
```

### 2. Run the simulation

Open `sim_cwxco.R` and run the script. By default it executes **1,000 iterations** with a sample size of **1,002 per group** using 8 parallel workers:

```r
plan(multisession, workers = 8)
niter <- 1000
sample_size <- 1002
```

Adjust these as needed for your machine. A single serial run for testing can be done with:

```r
sim_cwxco(N = 1001, prop_high_edu = .30, b0 = 0, b1 = .2, dem_prob_cut = .9)
```

### 3. Inspect results

`sim_cwxco()` returns a data frame comparing cocalibration and crosswalk regression coefficients across all scenarios. After parallel execution, results are collected as:

```r
results_df <- rbindlist(results, idcol = "rep") %>%
  mutate(n_sample = sample_size)
```

Cocalibration model objects are saved to `models/cc_models.Rds` on the first run and reused thereafter.

---

## Key Functions (`src/helper.R`)

| Function | Description |
|----------|-------------|
| `simCog()` | Simulates IRT-based cognitive item responses with a latent regression on education. |
| `cocalibrate()` | Places two groups on a common IRT scale using shared anchor items. |
| `getscore()` | Fits a unidimensional graded response IRT model and appends factor scores to a dataset. |
| `getcrosswalk()` | Wraps `cogxwalkr::crosswalk()` to map scores between two cognitive measures, conditioned on dementia status. |
| `cocalibrated_regressions()` | Fits education–cognition regressions on cocalibrated scores and extracts coefficients. |
| `get_scenario_data()` | Extracts scenario-specific item subsets from the full dataset. |
| `sim_cwxco()` | Master function: simulates data, runs cocalibration and crosswalks, and returns combined results. |

---

## Simulation Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `N` | 1002 | Sample size per group |
| `prop_high_edu` | 0.30 | Proportion with high education |
| `b0` | 0 | Intercept for the latent education regression |
| `b1` | 0.20 | Slope (true education effect on cognition) |
| `dem_prob_cut` | 0.90 | Probability threshold for assigning dementia status |
| `dem_on_edu` | `FALSE` | Whether dementia probability model includes education |
| `niter` | 1000 | Number of simulation replications |

---

## Notes

- Parallel execution uses the [`future`](https://future.futureverse.org/) framework. On HPC/cluster environments, replace `multisession` with an appropriate `future` plan (e.g., `future.batchtools`).
- Setting `future.globals.maxSize` may be necessary for large simulations (currently set to 1,500 MB).
- Memory usage scales with `N × niter`; reduce these for testing on lower-memory machines.
- Cocalibration models (`cc_models.Rds`) are only written once. Delete the file to regenerate them.
- A random seed is set per iteration (`set.seed(iter)`) to ensure reproducibility across runs.

