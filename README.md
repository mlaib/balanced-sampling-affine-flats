# Exactly balanced sampling on overlapping stratifications with affine flats

Code, data and results for the manuscript of the same title by Mohamed Laib. Every number in the
tables of the manuscript is produced by the scripts below from the files in `data/`, and every
intermediate result is stored in `results/`.

## Requirements

- Python 3.12 with numpy 2.4.6 and scipy 1.18.0 (`requirements.txt`), for all exact results.
- R 4.3.3 with BalancedSampling 2.1.1, for the Monte Carlo runs of the cube, local pivotal and local cube methods, and with sampling 2.11, for those of `samplecube`.
- A C compiler (gcc), for the enumeration of orthogonal arrays in `settle_attainment.py` and its cross-check in `crosscheck_oa16.py`.

`environment.txt` records the versions used for the stored results.

## Reproduction

```
bash reproduce.sh            # everything, including the Monte Carlo runs (several hours)
QUICK=1 bash reproduce.sh    # exact computations only, then the tables from the stored Monte Carlo results
```

`make_tables.py` is the entry point for the tables: it rebuilds `results/table1.md` to `results/table3.md`,
`results/tableA.1.md` and `results/tableS1.md` to `results/tableS7.md` (Tables 1 to 3, A.1 and S1 to S7 of the journal version of the manuscript) from the stored results. Release v1.9 wrote the numbering of the shortened version (Tables 2 to 5, A.1 and S1 to S7) and release v1.8 that of the version before it (Tables 2 to 7, A.1, S1 and S2).

## Scripts

| File | Purpose |
|---|---|
| `affsamp.py` | Linear algebra over GF(p), subspace enumeration, admissible directions, coset designs, spectral weights, exact variances, inclusion probabilities. |
| `test_theory.py` | Brute-force checks (support enumeration) of every exact statement in Sections 2–4, including Proposition 3; prints the number of checks per statement. |
| `study_factorial.py` | Section 5.1 (complete factorial frames 2^7 and 3^5): the five responses, exact design effects, minimax mixtures; writes the inputs of `mc_cube.R`. |
| `study_lattice.py` | Section 5.1: lattice sampling with a linear orthogonal array of maximal strength (Table 2, lattice column). |
| `mc_cube.R` | Section 5.1: Monte Carlo for the cube method (random unit order before each draw), with the split of each squared error into main-effect, residual and cross parts, the main-effect weight t_M, and optionally the selected units of every draw. |
| `ht_factorial.py` | Section 5.1, Tables 2, S3, S4 and S6: the Horvitz–Thompson mean (fixed denominator k) of the cube variants (BalancedSampling, its continued flight, both landings of `samplecube`, and the runs with the balancing variables in reverse order) from their stored samples: draws of size other than k, design effects of the Horvitz–Thompson and of the realised sample mean, the parts of the squared error due to the size, the main effects, the residual and the cross terms, and the main-effect weight share τ with its standard error. |
| `cube_factorial.py` | Section 5.1: the Python transcription of `cube` on the factorial frames, with the package's landing and with the flight continued to a vertex (20 000 draws per configuration). |
| `replicate_populations.py` | Supplementary material S3, Table S1: the factorial comparisons for 500 further populations (exact for coset designs and lattice sampling, stored samples for the cube variants). |
| `study_minimax_lp.py` | Table 1: minimax mixing distributions by linear programming, reduced by the monomial group. |
| `certify_lp.py` | Exact rational certificates for program (8): primal mixtures where the bound is attained, dual solutions where it is not. |
| `settle_attainment.py`, `c/enum_oa16.c` | Proposition 6 (Section 4): enumeration of all balanced samples where (8) stays above the bound; every one is a flat except for OA(16,7,2,2), which `c/enum_oa16.c` lists in a canonical form (`enum_oa16 t [m]`) and for which the minimax program over all arrays is solved and certified in rational arithmetic (compiles the C program with gcc). |
| `crosscheck_oa16.py`, `c/enum_direct.c` | Appendix A: independent check of that enumeration. `c/enum_direct.c` lists every balanced sample containing the origin by backtracking over the units, without symmetry reduction; its counts are compared, profile by profile, with those implied by the canonical forms, for m = 5, 6 and 7 (about half an hour). |
| `check_certificates.py` | Appendix A: standalone check, with the Python standard library only, of the certificates of Table 1 and of `results/attainment.json`; it lists the admissible directions afresh. |
| `study_grid.py` | Section 5.2 (61 x 61 elevation grid): exact design effects, ranking of the lattice directions, exact evaluation of the variance estimator (Table S7), residual Fourier energy; writes the inputs of `mc_grid.R`. |
| `mc_grid.R` | Section 5.2: Monte Carlo for the cube, local pivotal (lpm2, planar or wrapped distances) and local cube (lcube) methods, with the main-effect weights t_M on S2 and S4, the size and the Horvitz–Thompson estimate of every draw, and optionally the selected units. |
| `run_grid_mc.sh` | Runs `mc_grid.R` for every design and sample size of Section 5.2 with the seeds and replicate numbers of the supplementary material (S2). |
| `run_grid_mc_v2.sh` | Reruns with the same seeds for the main-effect weights (first 100 draws for the slow cube S4 and local cube runs) and the local pivotal method with wrapped distances; superseded by `run_grid_mc_v3.sh`. |
| `run_grid_mc_v3.sh` | Section 5.2, Tables 3 and S5: every grid design again with all its draws and the original seeds, recording the size, the Horvitz–Thompson estimate and the main-effect weights of every draw, and the selected units of the cube and local cube methods. |
| `cube_grid_fixed.py` | Section 5.2, Table 3: the cube method balanced on S4 with each direction computed from a QR factorisation (otherwise the transcription of `audit_cube.py`), 400 draws per sample size; it keeps the sample size, which the package loses on S4. |
| `check_rref_defect.py` | Supplementary material S4: the transcription with the package's row reduction (exact zero tests) and with a tolerant one, from the same seeds: every factorial draw whose stored size differed from k, 300 control draws per configuration, and a 19 x 19 grid with two to four line stratifications. |
| `check_size_drift.R` | Supplementary material S4: sizes of the samples of the package's `cube` and `lcube` on a 19 x 19 grid balanced on two to four line stratifications. |
| `bound_centred.py` | Section 5 and Tables 3, S2 and S5: the lower bound (S1) of the supplement, Proposition 3(a) for the Monte Carlo designs, corrected for a varying sample size (responses with mean zero), from the per-draw sizes and stratum counts, with a joint standard error. |
| `mc_bias.py` | Section 5 and Section S2: the mean error of the Horvitz–Thompson estimator, with its Monte Carlo standard error, for every stored Monte Carlo sample (252 combinations of design, sample size and response), including the draws of another size. |
| `check_ranks.py` | Supplementary material S6: the ranks of the sample main-effects model used by the regression variance estimator, by singular values and exactly by elimination modulo a prime, for every sample of the factorial supports. |
| `bootstrap_bank.py` | Supplementary material S3: resamples whole draws from the stored bank of each cube variant, the same resample for every population, to give the Monte Carlo standard errors of the cube quantiles of Table S1. |
| `varest_exact.py` | Section 4 and Supplement S6: exact evaluation of the random-group and regression variance estimators (factorial replicated designs, grid design balanced on rows). |
| `audit_cube.py` | Section 5.2: Python transcription of the flight and landing phases of `cube` in BalancedSampling 2.1.1, instrumented, with the flight optionally continued until no move preserves all balancing equations; `--decompose` splits each mean squared error into main-effect, residual and cross parts. Given a sequence of uniform numbers (`UniformStream`), it uses them where the package calls R's generator. |
| `validate_transcription.py` | Supplementary material S4: runs `cube` in R for ten seeds per sample size (grid and factorial frames), feeds the uniform numbers R generates for each seed to `audit_cube.py`, and checks that both return the same sample. |
| `count_balanced_samples.py` | Sections 2 and 4 and Supplementary material S1: vertices and integral vertices of the balancing polytope of Example 2 of Tillé and Rivest (N = 8 and 16), and the samples of the q x q grid balanced on rows, columns and anti-diagonals (q = 5, 7, 11), with the weights of the uniform design over them. |
| `time_preprocessing.py` | Section 4: time and memory to enumerate admissible directions and solve program (8), acceptance rates of rejection sampling, size of the stored minimax supports. |
| `make_tables.py` | Tables 1–3, A.1 and S1–S7 from `results/`, including the exact shares of grid coset samples balanced on S2 and S4 and the exact main-effect weights of the coset designs. |
| `run_samplecube.sh` | Section 5 and Supplement S2: `samplecube` of the package sampling 2.11 through `mc_cube.R` and `mc_grid.R` (both landings in the factorial frames, the landing by suppression of variables on S2 and S4 of the grid), and BalancedSampling and both landings of `samplecube` with the balancing variables in reverse order; each job is skipped if its output exists, so an interrupted run can be resumed. |
| `samplecube_grid.py` | Section 5.2, Tables 3 and S5: summaries of the grid draws of `samplecube`, combined over chunks: sizes, design effects, shares of exactly balanced samples and τ. |
| `multisurface/` | Supplementary material S8: the best-L coset designs of Section 5.2 on 256 further 61 x 61 elevation blocks from 16 terrains and on 360 simulated Gaussian fields, against the local pivotal, local cube and cube methods. `PLAN.md` holds the designs, surfaces and decision rule, written before any result was computed, and the outcome; `scripts/fetch_tiles.py` downloads the terrain tiles and cuts the blocks, `scripts/simulate_fields.py` simulates the fields, `scripts/draw_mc.R` and `scripts/run_mc.sh` store the Monte Carlo samples (evaluated on every surface), and `scripts/evaluate.py` writes `results/deff_long.csv`, `results/summary.md` and `results/tableS8.md`. `scripts/explore_systematic_cosets.py` is an exploratory check made after the results (systematic choice of the cosets for k = 183) and is not used in the manuscript. |
| `fig_coset_samples.py` | Figure 1 (Section 5.2): coset samples of the last- and first-ranked lattice directions on the 61 x 61 grid, over elevation contours; writes `results/fig1_coset_samples.pdf`, `.png` and a compact copy for BrainLog. |
| `test_prime_power.py` | Supplement S1: a prime power reduces to its prime field (GF(4)^2 as GF(2)^4 with two and three GF(4)-line stratifications); writes `results/prime_power_check.txt` when its output is redirected there. |
| `reproduce.sh` | Runs all of the above in order. |

## Data

| File | Content |
|---|---|
| `multisurface/data/elevation_blocks.npz` | The 256 elevation blocks of Supplement S8: Terrain Tiles on the Registry of Open Data on AWS (Mapzen/Tilezen, Terrarium encoding, zoom 10), with sources including SRTM, 3DEP and GMTED2010 data courtesy of the U.S. Geological Survey and EU-DEM layers produced using Copernicus data and information funded by the European Union; see https://github.com/tilezen/joerd/blob/master/docs/attribution.md for the full list and licences. |
| `data/volcano61.csv` | Rows 14–74 (all 61 columns) of R's `datasets::volcano`; unit (i, j) has code 61 i + j. |
| `data/grid_y.csv`, `data/grid_X.csv` | Elevations in unit order, and the 0/1 indicators of rows, columns, anti-diagonals i+j and diagonals i−j (levels 1–60 of each; level 0 dropped). |
| `data/fact_p2_m7_Y.csv`, `data/fact_p3_m5_Y.csv` | The five responses of Section 5.1 (additive, mixed, interaction, noise, higher-order), seeds 7 and 5. |
| `data/fact_p2_m7_X.csv`, `data/fact_p3_m5_X.csv` | Category indicators used as balancing variables by the cube method. |

## Results

Directions are stored as reduced row-echelon bases of the dual subspace D = V⊥ of the flat direction V.

| File | Content |
|---|---|
| `results/factorial_exact.csv` | Exact design effects of the coset designs of Section 5.1, with R² and adjusted R². |
| `results/factorial_minimax_mixtures.json` | Support (dual bases) and probabilities of the minimax mixtures of Section 5.1. |
| `results/factorial_lattice.csv` | Lattice sampling: strength, dual basis of the chosen array, exact design effects and worst case. |
| `results/factorial_designs.json` | Size of the support, worst case relative to the bound, positivity of second-order inclusion probabilities. |
| `results/factorial_cube_p2_m7.csv`, `results/factorial_cube_p3_m5.csv` | Monte Carlo summaries for the cube method (mean squared error, standard error, share exactly balanced, the main-effect, residual and cross parts of the mean squared error with their standard errors, and the main-effect weight). |
| `results/cube_samples/factorial_*.txt.gz`, `results/cube_samples/transcription_*.txt.gz` | Selected units of every factorial draw: package, and transcription with both landings. |
| `results/factorial_ht.csv` | Horvitz–Thompson results of the factorial cube variants: draws of size other than k, design effects of the Horvitz–Thompson and sample means, parts of the squared error, shares exactly balanced and τ, with standard errors. |
| `results/cube_factorial.csv` | Transcription on the factorial frames: design effects, shares exactly balanced, shares ending at an integral vertex, main-effect weights. |
| `results/replicated_populations.csv`, `results/replicated_populations.npz` | Quantiles of design effects and of ratios to the minimax design over 500 populations, and all design effects. |
| `results/attainment.json` | Enumeration of balanced samples where (8) stays above the bound; profiles and exact linear program for OA(16,7,2,2). |
| `results/crosscheck_oa16.json` | Counts of the canonical and direct enumerations of OA(16,m,2,2) for m = 5, 6, 7 and the outcome of each check. |
| `results/varest_exact.csv` | Exact moments and coverage of the variance estimators. |
| `results/minimax_lp.csv` | Value of program (8) for every configuration solved (Table 1 shows a subset). |
| `results/lp_certificates.json` | Rational certificates: mixtures attaining the bound, and dual solutions with the exact minimax value where it is not attained. |
| `results/grid_exact.json` | Exact design effects of the grid designs, the evaluation of the variance estimator, and the sensitivity to ties among lattice directions. |
| `results/grid_lattice_ranking.csv` | The 58 directions admissible for S4 ranked by the length of the shortest dual vector; ties broken by increasing slope. |
| `results/grid_spectrum.json` | Residual Fourier energy by wrapped frequency length and on the dual lines of the chosen designs; seams of the wrapped grid. |
| `results/grid_balance_shares.json` | Exact share of the samples of each grid coset design that are balanced on S2 and on S4. |
| `results/mc_grid/*.csv` | Monte Carlo summaries for the cube, local pivotal and local cube methods on the grid. |
| `results/mc_grid_v2/*.csv` | Reruns with main-effect weights, and the local pivotal method with wrapped distances. |
| `results/mc_grid_v3/*.csv` | Reruns of every grid design with all draws (Tables 3 and S5): mean squared errors of the sample mean and of the Horvitz–Thompson mean, draws of size other than k, τ on S2 and S4 with standard errors. |
| `results/cube_samples/grid_*.txt.gz` | Selected units of every grid draw of the cube and local cube methods (one line per draw: k, design, units). |
| `results/cube_grid_fixed.csv`, `results/cube_grid_fixed_draws.csv` | The cube method on S4 with directions by QR factorisation: summary (design effects, sizes, shares balanced, τ) and one line per draw. |
| `results/rref_defect_check.csv`, `results/size_drift_check.csv` | Sizes returned with the package's and with a tolerant row reduction (transcription), and by the package on a 19 x 19 grid. |
| `results/cube_audit_S2.csv`, `results/cube_audit_S2_draws.csv` | Audit of the cube method balanced on rows and columns: summary, and one line per draw with the selected units. |
| `results/cube_audit_S2_decomposition.csv` | Decomposition of the mean squared error of the audited designs into main-effect, residual and cross parts. |
| `results/exhaustive_counts.json` | Exhaustive counts of Sections 1 and 4: vertices of the Example 2 polytopes, balanced samples of the q x q grids, flats among them and spectral weights of the uniform design. |
| `results/transcription_check.csv` | For each sample size and seed, whether the transcription returned the package's sample, and how many uniform numbers it used. |
| `results/preprocessing_costs.csv` | Preprocessing cost of the coset designs of Table 1 and Section 5.1 (one laptop core). |
| `results/bound_centred.csv` | The lower bounds (S1) of the Monte Carlo designs, with the size correction, the main-effect weight share and their standard errors. |
| `results/mc_bias.csv` | The mean errors of the Horvitz–Thompson estimator over the stored samples, with standard errors, z values and the number of draws of another size. |
| `results/rank_check.json` | Singular-value gap and exact ranks of the sample design matrices of the regression estimator. |
| `results/bootstrap_bank.csv` | Bootstrap standard errors of the cube quantiles and shares of Table S1, and the Monte Carlo error of the ratio within a population. |
| `results/factorial_scube1_*.csv`, `results/factorial_scube2_*.csv`, `results/factorial_*_rev_*.csv` | Monte Carlo summaries of `samplecube` (landing by linear programming and by suppression of variables) and of the runs with the balancing variables in reverse order; their selected units are in `results/cube_samples/`. |
| `results/mc_grid_scube/*.csv`, `results/samplecube_grid.csv` | Monte Carlo summaries of `samplecube` on the grid, by chunk and combined. |
| `results/table1.md` – `results/table3.md`, `results/tableA.1.md`, `results/tableS1.md` – `results/tableS7.md` | Tables 1 to 3, A.1 and S1 to S7 of the manuscript (the file name is the table number). |
