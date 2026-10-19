# Exactly balanced sampling on overlapping stratifications with affine flats: code and data

[![License: AGPL-3.0](https://img.shields.io/badge/license-AGPL--3.0-blue.svg)](LICENSE)
[![Python 3.12](https://img.shields.io/badge/Python-3.12-3776AB.svg)](requirements.txt)
[![R 4.3](https://img.shields.io/badge/R-4.3-276DC3.svg)](environment.txt)

This repository contains the code, data and stored results of the manuscript *Exactly balanced sampling on
overlapping stratifications with affine flats* by Mohamed Laib. Every number in the tables and figures of the
manuscript and its supplementary material is produced by the scripts here. Every intermediate result, including
the selected units of every Monte Carlo draw, is stored in `results/`, so the tables can be rebuilt without
rerunning the simulations.

The designs are also available as an R package, [cosetsampling](https://github.com/mlaib/cosetsampling).

## Contents

```
affsamp.py          core library: linear algebra over GF(p), admissible directions, coset designs,
                    spectral weights, exact variances and inclusion probabilities
*.py, *.R, *.sh     one script per analysis (see "From the manuscript to the code")
c/                  C programs that enumerate the orthogonal arrays OA(16, m, 2, 2)
data/               populations, responses and balancing variables
multisurface/       Supplement S8: further elevation blocks and simulated fields
results/            stored results, Monte Carlo samples and the tables of the manuscript
reproduce.sh        runs every analysis in order
```

## Requirements

- **Python 3.12** with numpy 2.4.6 and scipy 1.18.0 for all exact results, and matplotlib and Pillow for
  Figure 1 and the download of the terrain tiles: `pip install -r requirements.txt`.
- **R 4.3.3** with BalancedSampling 2.1.1 and sampling 2.11, for the Monte Carlo runs only. The stored results
  depend on these versions, for example with
  `remotes::install_version("BalancedSampling", "2.1.1")` and `remotes::install_version("sampling", "2.11")`.
- **A C compiler** (gcc) for the enumerations in `c/`, which `settle_attainment.py` and `crosscheck_oa16.py`
  compile and run.

`environment.txt` records the software versions used for the stored results.

## Reproducing the results

```bash
pip install -r requirements.txt
QUICK=1 bash reproduce.sh   # exact computations and enumerations, then all tables from the stored Monte Carlo samples
bash reproduce.sh           # everything, including the Monte Carlo simulations (several hours)
```

The Monte Carlo simulations take several hours; among the exact computations, the independent cross-check of the
OA(16,7,2,2) enumeration alone takes about half an hour. Every random computation has a fixed seed, recorded in
the script and in Section S2 of the supplementary material.

The tables are written to `results/`:
- `make_tables.py` writes Tables 1 to 3, A.1 and S1 to S7 (`results/table1.md`, …, `results/tableS7.md`);
- `multisurface/scripts/evaluate.py` writes Table S8 (`multisurface/results/tableS8.md`);
- `table_s9.py` writes Table S9 (`results/tableS9.md`);
- `fig_coset_samples.py` draws Figure 1 (`results/fig1_coset_samples.pdf`).

## From the manuscript to the code

| Manuscript | Scripts | Main outputs in `results/` |
|---|---|---|
| Sections 2 to 4: every exact statement | `test_theory.py`, `test_prime_power.py` | checks printed; `prime_power_check.txt` |
| Sections 2 and 4, Supplement S1: exhaustive counts | `count_balanced_samples.py` | `exhaustive_counts.json` |
| Table 1: minimax mixtures and their certificates | `study_minimax_lp.py`, `certify_lp.py` | `minimax_lp.csv`, `lp_certificates.json`, `table1.md` |
| Proposition 6, Appendix A, Tables A.1 and S9 | `settle_attainment.py`, `c/enum_oa16.c`, `crosscheck_oa16.py`, `c/enum_direct.c`, `check_certificates.py`, `table_s9.py` | `attainment.json`, `crosscheck_oa16.json`, `tableA.1.md`, `tableS9.md` |
| Section 4: preprocessing cost and variance estimators; Supplement S6, Table S7 | `time_preprocessing.py`, `varest_exact.py`, `check_ranks.py` | `preprocessing_costs.csv`, `varest_exact.csv`, `rank_check.json`, `tableS7.md` |
| Section 5.1: complete factorial frames; Tables 2, S3, S4 and S6 | `study_factorial.py`, `study_lattice.py`, `mc_cube.R`, `cube_factorial.py`, `run_samplecube.sh`, `ht_factorial.py` | `factorial_*.csv`, `factorial_*.json`, `cube_samples/`, `table2.md` |
| Supplement S3, Table S1: 500 further populations | `replicate_populations.py`, `bootstrap_bank.py` | `replicated_populations.*`, `bootstrap_bank.csv`, `tableS1.md` |
| Section 5.2: the 61 × 61 elevation grid; Tables 3, S2 and S5; Figure 1 | `study_grid.py`, `mc_grid.R`, `run_grid_mc_v3.sh`, `cube_grid_fixed.py`, `samplecube_grid.py`, `bound_centred.py`, `fig_coset_samples.py` | `grid_*.json`, `grid_lattice_ranking.csv`, `mc_grid_v3/`, `cube_grid_fixed.csv`, `table3.md`, `fig1_coset_samples.pdf` |
| Section 5, Supplement S2: bias of the estimator and centring of the response | `mc_bias.py` | `mc_bias.csv` |
| Supplement S4: diagnostics of the implementation of the cube method | `audit_cube.py`, `validate_transcription.py`, `check_rref_defect.py`, `check_size_drift.R` | `cube_audit_S2*.csv`, `transcription_check.csv`, `rref_defect_check.csv`, `size_drift_check.csv` |
| Section 5.3, Supplement S8, Table S8: further elevation blocks and simulated fields | `multisurface/` | `multisurface/results/` |

<details>
<summary><b>Description of every script</b></summary>

| File | Purpose |
|---|---|
| `affsamp.py` | Linear algebra over GF(p), subspace enumeration, admissible directions, coset designs, spectral weights, exact variances, inclusion probabilities. |
| `test_theory.py` | Brute-force checks, by enumeration of the support, of every exact statement in Sections 2 to 4, including Proposition 3; prints the number of checks per statement. |
| `test_prime_power.py` | Supplement S1: a prime power reduces to its prime field (GF(4)^2 as GF(2)^4 with two and three GF(4)-line stratifications); its output is stored in `results/prime_power_check.txt`. |
| `count_balanced_samples.py` | Sections 2 and 4, Supplement S1: vertices and integral vertices of the balancing polytope of Example 2 of Tillé and Rivest (N = 8 and 16), and the samples of the q × q grid balanced on rows, columns and anti-diagonals (q = 5, 7, 11), with the weights of the uniform design over them. |
| `study_minimax_lp.py` | Table 1: minimax mixing distributions by linear programming, reduced by the monomial group. |
| `certify_lp.py` | Exact rational certificates for program (8): primal mixtures where the bound is attained, dual solutions where it is not. |
| `settle_attainment.py`, `c/enum_oa16.c` | Proposition 6: enumeration of all balanced samples where (8) stays above the bound. Every one is a flat except for OA(16,7,2,2), which `c/enum_oa16.c` lists in a canonical form (`enum_oa16 t [m]`) and for which the minimax program over all arrays is solved and certified in rational arithmetic. |
| `crosscheck_oa16.py`, `c/enum_direct.c` | Appendix A: independent check of that enumeration. `c/enum_direct.c` lists every balanced sample containing the origin by backtracking over the units, without symmetry reduction; its counts are compared, profile by profile, with those implied by the canonical forms, for m = 5, 6 and 7 (about half an hour). |
| `check_certificates.py` | Appendix A and Supplement S9: standalone check of the certificates of Table 1 and Appendix A with the Python standard library only; it lists the admissible directions afresh. |
| `table_s9.py` | Supplement S9, Table S9: the verification map of Proposition 6, from the stored certificates and enumerations. |
| `time_preprocessing.py` | Section 4: time and memory to enumerate admissible directions and solve program (8), acceptance rates of rejection sampling, size of the stored minimax supports. |
| `varest_exact.py` | Section 4 and Supplement S6: exact evaluation of the random-group and regression variance estimators (factorial replicated designs, grid design balanced on rows). |
| `check_ranks.py` | Supplement S6: the ranks of the sample main-effects model used by the regression variance estimator, by singular values and exactly by elimination modulo a prime, for every sample of the factorial supports. |
| `study_factorial.py` | Section 5.1 (complete factorial frames 2^7 and 3^5): the five responses, exact design effects, minimax mixtures; writes the inputs of `mc_cube.R`. |
| `study_lattice.py` | Section 5.1: lattice sampling with a linear orthogonal array of maximal strength (Table 2, lattice column). |
| `mc_cube.R` | Section 5.1: Monte Carlo for the cube method (random unit order before each draw), with the split of each squared error into main-effect, residual and cross parts, the main-effect weight t_M, and optionally the selected units of every draw. |
| `cube_factorial.py` | Section 5.1: the Python transcription of `cube` on the factorial frames, with the package's landing and with the flight continued to a vertex (20 000 draws per configuration). |
| `run_samplecube.sh` | Section 5 and Supplement S2: `samplecube` of the package sampling 2.11 through `mc_cube.R` and `mc_grid.R` (both landings in the factorial frames, the landing by suppression of variables on S2 and S4 of the grid), and BalancedSampling and both landings of `samplecube` with the balancing variables in reverse order. Each job is skipped if its output exists, so an interrupted run can be resumed. |
| `ht_factorial.py` | Tables 2, S3, S4 and S6: the Horvitz–Thompson mean (fixed denominator k) of the cube variants from their stored samples: draws of size other than k, design effects of the Horvitz–Thompson and of the realised sample mean, the parts of the squared error due to the size, the main effects, the residual and the cross terms, and the main-effect weight share τ with its standard error. |
| `replicate_populations.py` | Supplement S3, Table S1: the factorial comparisons for 500 further populations (exact for coset designs and lattice sampling, stored samples for the cube variants). |
| `bootstrap_bank.py` | Supplement S3: resamples whole draws from the stored bank of each cube variant, the same resample for every population, to give the Monte Carlo standard errors of the cube quantiles of Table S1. |
| `study_grid.py` | Section 5.2 (61 × 61 elevation grid): exact design effects, ranking of the lattice directions, exact evaluation of the variance estimator (Table S7), residual Fourier energy; writes the inputs of `mc_grid.R`. |
| `mc_grid.R` | Section 5.2: Monte Carlo for the cube, local pivotal (`lpm2`, planar or wrapped distances) and local cube (`lcube`) methods, with the main-effect weights t_M on S2 and S4, the size and the Horvitz–Thompson estimate of every draw, and optionally the selected units. |
| `run_grid_mc_v3.sh` | Section 5.2, Tables 3 and S5: every grid design with all its draws and the original seeds, recording the size, the Horvitz–Thompson estimate and the main-effect weights of every draw, and the selected units of the cube and local cube methods. The earlier runs with the same seeds, `run_grid_mc.sh` and `run_grid_mc_v2.sh`, are kept for the record. |
| `cube_grid_fixed.py` | Section 5.2, Table 3: the cube method balanced on S4 with each direction computed from a QR factorisation (otherwise the transcription of `audit_cube.py`), 400 draws per sample size; it keeps the sample size, which the package loses on S4. |
| `samplecube_grid.py` | Section 5.2, Tables 3 and S5: summaries of the grid draws of `samplecube`, combined over chunks: sizes, design effects, shares of exactly balanced samples and τ. |
| `bound_centred.py` | Section 5, Tables 3, S2 and S5: the lower bound (S1) of the supplement, Proposition 3(a) for the Monte Carlo designs, corrected for a varying sample size (responses with mean zero), from the per-draw sizes and stratum counts, with a joint standard error. |
| `mc_bias.py` | Section 5 and Supplement S2: the mean error of the Horvitz–Thompson estimator, with its Monte Carlo standard error, for every stored Monte Carlo sample (252 combinations of design, sample size and response), including the draws of another size, and the design effects for the response and for the centred response. |
| `fig_coset_samples.py` | Figure 1 (Section 5.2): coset samples of the last- and first-ranked lattice directions on the 61 × 61 grid, over elevation contours; writes `results/fig1_coset_samples.pdf` and `.png`. |
| `audit_cube.py` | Supplement S4: Python transcription of the flight and landing phases of `cube` in BalancedSampling 2.1.1, instrumented, with the flight optionally continued until no move preserves all balancing equations; `--decompose` splits each mean squared error into main-effect, residual and cross parts. Given a sequence of uniform numbers (`UniformStream`), it uses them where the package calls R's generator. |
| `validate_transcription.py` | Supplement S4: runs `cube` in R for ten seeds per sample size (grid and factorial frames), feeds the uniform numbers R generates for each seed to `audit_cube.py`, and checks that both return the same sample. |
| `check_rref_defect.py` | Supplement S4: the transcription with the package's row reduction (exact zero tests) and with a tolerant one, from the same seeds: every factorial draw whose stored size differed from k, 300 control draws per configuration, and a 19 × 19 grid with two to four line stratifications. |
| `check_size_drift.R` | Supplement S4: sizes of the samples of the package's `cube` and `lcube` on a 19 × 19 grid balanced on two to four line stratifications. |
| `multisurface/` | Supplement S8: the best-L coset designs of Section 5.2 on 256 further 61 × 61 elevation blocks from 16 terrains and on 360 simulated Gaussian fields, against the local pivotal, local cube and cube methods. `PLAN.md` holds the designs, surfaces and decision rule, written before any result was computed, and the outcome. `scripts/fetch_tiles.py` downloads the terrain tiles and cuts the blocks, `scripts/simulate_fields.py` simulates the fields, `scripts/draw_mc.R` and `scripts/run_mc.sh` store the Monte Carlo samples, which are evaluated on every surface, and `scripts/evaluate.py` writes `results/deff_long.csv`, `results/summary.md` and `results/tableS8.md`. `scripts/explore_systematic_cosets.py` is an exploratory check made after the results and is not used in the manuscript. |
| `make_tables.py` | Tables 1 to 3, A.1 and S1 to S7 from `results/`, including the exact shares of grid coset samples balanced on S2 and S4 and the exact main-effect weights of the coset designs. |
| `reproduce.sh` | Runs all of the above in order. |

</details>

## Data

| File | Content |
|---|---|
| `data/volcano61.csv` | Rows 14 to 74 (all 61 columns) of R's data set `datasets::volcano`, the topographic map of Maunga Whau; unit (i, j) has code 61 i + j. |
| `data/grid_y.csv`, `data/grid_X.csv` | Elevations in unit order, and the 0/1 indicators of rows, columns, anti-diagonals i + j and diagonals i − j (levels 1 to 60 of each; level 0 dropped). |
| `data/fact_p2_m7_Y.csv`, `data/fact_p3_m5_Y.csv` | The five responses of Section 5.1 (additive, mixed, interaction, noise, higher-order), generated with seeds 7 and 5. |
| `data/fact_p2_m7_X.csv`, `data/fact_p3_m5_X.csv` | Category indicators used as balancing variables by the cube method. |
| `multisurface/data/elevation_blocks.npz` | The 256 elevation blocks of Supplement S8, from the [Terrain Tiles](https://registry.opendata.aws/terrain-tiles/) on the Registry of Open Data on AWS (Mapzen/Tilezen, Terrarium encoding, zoom 10). |

The Terrain Tiles combine sources that include SRTM, 3DEP and GMTED2010 data courtesy of the U.S. Geological
Survey, and EU-DEM layers produced using Copernicus data and information funded by the European Union. The full
list of sources and their licences is in the
[Tilezen attribution](https://github.com/tilezen/joerd/blob/master/docs/attribution.md).

## Stored results

Directions are stored as reduced row-echelon bases of the dual subspace D = V⊥ of the flat direction V.

<details>
<summary><b>Description of every result file</b></summary>

| File | Content |
|---|---|
| `results/table1.md` to `results/table3.md`, `results/tableA.1.md`, `results/tableS1.md` to `results/tableS7.md`, `results/tableS9.md` | Tables 1 to 3, A.1, S1 to S7 and S9 of the manuscript (the file name is the table number); Table S8 is `multisurface/results/tableS8.md`. |
| `results/factorial_exact.csv` | Exact design effects of the coset designs of Section 5.1, with R² and adjusted R². |
| `results/factorial_minimax_mixtures.json` | Support (dual bases) and probabilities of the minimax mixtures of Section 5.1. |
| `results/factorial_lattice.csv` | Lattice sampling: strength, dual basis of the chosen array, exact design effects and worst case. |
| `results/factorial_designs.json` | Size of the support, worst case relative to the bound, positivity of second-order inclusion probabilities. |
| `results/factorial_cube_p2_m7.csv`, `results/factorial_cube_p3_m5.csv` | Monte Carlo summaries for the cube method: mean squared error and its standard error, share exactly balanced, the main-effect, residual and cross parts of the mean squared error with their standard errors, and the main-effect weight. |
| `results/factorial_scube1_*.csv`, `results/factorial_scube2_*.csv`, `results/factorial_*_rev_*.csv` | Monte Carlo summaries of `samplecube` (landing by linear programming and by suppression of variables) and of the runs with the balancing variables in reverse order. |
| `results/cube_samples/` | Selected units of every Monte Carlo draw of the cube variants, the local cube method and `samplecube`, in the factorial frames and on the grid (one line per draw). |
| `results/factorial_ht.csv` | Horvitz–Thompson results of the factorial cube variants: draws of size other than k, design effects of the Horvitz–Thompson and sample means, parts of the squared error, shares exactly balanced and τ, with standard errors. |
| `results/cube_factorial.csv` | The transcription on the factorial frames: design effects, shares exactly balanced, shares ending at an integral vertex, main-effect weights. |
| `results/replicated_populations.csv`, `results/replicated_populations.npz` | Quantiles of design effects and of ratios to the minimax design over 500 populations, and all design effects. |
| `results/bootstrap_bank.csv` | Bootstrap standard errors of the cube quantiles and shares of Table S1, and the Monte Carlo error of the ratio within a population. |
| `results/minimax_lp.csv` | Value of program (8) for every configuration solved (Table 1 shows a subset). |
| `results/lp_certificates.json` | Rational certificates: mixtures attaining the bound, and dual solutions with the exact minimax value where it is not attained. |
| `results/attainment.json` | Enumeration of balanced samples where (8) stays above the bound; profiles and exact linear program for OA(16,7,2,2). |
| `results/crosscheck_oa16.json` | Counts of the canonical and direct enumerations of OA(16, m, 2, 2) for m = 5, 6, 7 and the outcome of each check. |
| `results/exhaustive_counts.json` | Exhaustive counts of Sections 2 and 4: vertices of the Example 2 polytopes, balanced samples of the q × q grids, flats among them and spectral weights of the uniform design. |
| `results/prime_power_check.txt` | Output of `test_prime_power.py`. |
| `results/varest_exact.csv` | Exact moments and coverage of the variance estimators. |
| `results/rank_check.json` | Singular-value gap and exact ranks of the sample design matrices of the regression estimator. |
| `results/preprocessing_costs.csv` | Preprocessing cost of the coset designs of Table 1 and Section 5.1, on one laptop core. |
| `results/grid_exact.json` | Exact design effects of the grid designs, the evaluation of the variance estimator, and the sensitivity to ties among lattice directions. |
| `results/grid_lattice_ranking.csv` | The 58 directions admissible for S4, ranked by the length of the shortest dual vector; ties broken by increasing slope. |
| `results/grid_spectrum.json` | Residual Fourier energy by wrapped frequency length and on the dual lines of the chosen designs; seams of the wrapped grid. |
| `results/grid_balance_shares.json` | Exact share of the samples of each grid coset design that are balanced on S2 and on S4. |
| `results/mc_grid_v3/*.csv` | Monte Carlo summaries of every grid design with all draws (Tables 3 and S5): mean squared errors of the sample mean and of the Horvitz–Thompson mean, draws of size other than k, τ on S2 and S4 with standard errors. The earlier runs are in `results/mc_grid/` and `results/mc_grid_v2/`. |
| `results/mc_grid_scube/*.csv`, `results/samplecube_grid.csv` | Monte Carlo summaries of `samplecube` on the grid, by chunk and combined. |
| `results/cube_grid_fixed.csv`, `results/cube_grid_fixed_draws.csv` | The cube method on S4 with directions by QR factorisation: summary (design effects, sizes, shares balanced, τ) and one line per draw. |
| `results/bound_centred.csv` | The lower bounds (S1) of the Monte Carlo designs, with the size correction, the main-effect weight share and their standard errors. |
| `results/mc_bias.csv` | Mean errors of the Horvitz–Thompson estimator over the stored samples, with standard errors, z values, the number of draws of another size, and the design effects for the response and for the centred response. |
| `results/cube_audit_S2.csv`, `results/cube_audit_S2_draws.csv`, `results/cube_audit_S2_decomposition.csv` | Audit of the cube method balanced on rows and columns: summary, one line per draw with the selected units, and the decomposition of the mean squared error into main-effect, residual and cross parts. |
| `results/transcription_check.csv` | For each sample size and seed, whether the transcription returned the package's sample, and how many uniform numbers it used. |
| `results/rref_defect_check.csv`, `results/size_drift_check.csv` | Sizes returned with the package's and with a tolerant row reduction (transcription), and by the package on a 19 × 19 grid. |
| `results/fig1_coset_samples.pdf`, `results/fig1_coset_samples.png` | Figure 1. |

</details>

## Citation

If you use this code or these data, please cite the manuscript:

> Laib, M. (2026). *Exactly balanced sampling on overlapping stratifications with affine flats*. Manuscript.

## Licence

GNU Affero General Public License v3.0 (AGPL-3.0). See [LICENSE](LICENSE).

## Contact

Mohamed Laib, Luxembourg Institute of Science and Technology (LIST).
Contact: laib.med@gmail.com, mohamed.laib@list.lu.
