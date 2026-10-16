# Exactly balanced sampling on overlapping stratifications with affine flats

Code, data and results for the manuscript of the same title by Mohamed Laib. Every number in the
tables of the manuscript is produced by the scripts below from the files in `data/`, and every
intermediate result is stored in `results/`.

## Requirements

- Python 3.12 with numpy 2.4.6 and scipy 1.18.0 (`requirements.txt`), for all exact results.
- R 4.3.3 with BalancedSampling 2.1.1, for the Monte Carlo runs of the cube, local pivotal and local cube methods.
- A C compiler (gcc), for the enumeration of orthogonal arrays in `settle_attainment.py`.

`environment.txt` records the versions used for the stored results.

## Reproduction

```
bash reproduce.sh            # everything, including the Monte Carlo runs (several hours)
QUICK=1 bash reproduce.sh    # exact computations only, then the tables from the stored Monte Carlo results
```

`make_tables.py` is the entry point for the tables: it rebuilds `results/table1.md` to `results/table4.md`
from the stored results.

## Scripts

| File | Purpose |
|---|---|
| `affsamp.py` | Linear algebra over GF(p), subspace enumeration, admissible directions, coset designs, spectral weights, exact variances, inclusion probabilities. |
| `test_theory.py` | Brute-force checks (support enumeration) of every exact statement in Sections 2–6, including Proposition 3; prints the number of checks per statement. |
| `study_factorial.py` | Section 7.1 (complete factorial frames 2^7 and 3^5): the five responses, exact design effects, minimax mixtures; writes the inputs of `mc_cube.R`. |
| `study_lattice.py` | Section 7.1: lattice sampling with a linear orthogonal array of maximal strength (Table 2, lattice column). |
| `mc_cube.R` | Section 7.1: Monte Carlo for the cube method (random unit order before each draw), with the split of each squared error into main-effect, residual and cross parts, the main-effect weight t_M, and optionally the selected units of every draw. |
| `cube_factorial.py` | Section 7.1: the Python transcription of `cube` on the factorial frames, with the package's landing and with the flight continued to a vertex (20 000 draws per configuration). |
| `replicate_populations.py` | Appendix A, Table A1: the factorial comparisons for 500 further populations (exact for coset designs and lattice sampling, stored samples for the cube variants). |
| `study_minimax_lp.py` | Table 1: minimax mixing distributions by linear programming, reduced by the monomial group. |
| `certify_lp.py` | Exact rational certificates for program (8): primal mixtures where the bound is attained, dual solutions where it is not. |
| `settle_attainment.py`, `c/enum_oa16.c` | Section 6: enumeration of all balanced samples where (8) stays above the bound; every one is a flat except for OA(16,7,2,2), for which the minimax program over all arrays is solved and certified in rational arithmetic (compiles the C program with gcc). |
| `study_grid.py` | Section 7.2 (61 x 61 elevation grid): exact design effects, ranking of the lattice directions, exact evaluation of the variance estimator (Table 4), residual Fourier energy; writes the inputs of `mc_grid.R`. |
| `mc_grid.R` | Section 7.2: Monte Carlo for the cube, local pivotal (lpm2, planar or wrapped distances) and local cube (lcube) methods, with the main-effect weights t_M on S2 and S4. |
| `run_grid_mc.sh` | Runs `mc_grid.R` for every design and sample size of Section 7.2 with the seeds and replicate numbers of Appendix A. |
| `run_grid_mc_v2.sh` | Reruns with the same seeds for the main-effect weights of Table 4 (first 100 draws for the slow cube S4 and local cube runs) and the local pivotal method with wrapped distances. |
| `varest_exact.py` | Section 7.3: exact evaluation of the random-group and regression variance estimators (factorial replicated designs, grid design balanced on rows). |
| `audit_cube.py` | Section 7.2: Python transcription of the flight and landing phases of `cube` in BalancedSampling 2.1.1, instrumented, with the flight optionally continued until no move preserves all balancing equations; `--decompose` splits each mean squared error into main-effect, residual and cross parts. Given a sequence of uniform numbers (`UniformStream`), it uses them where the package calls R's generator. |
| `validate_transcription.py` | Appendix A: runs `cube` in R for ten seeds per sample size (grid and factorial frames), feeds the uniform numbers R generates for each seed to `audit_cube.py`, and checks that both return the same sample. |
| `count_balanced_samples.py` | Sections 1 and 6: vertices and integral vertices of the balancing polytope of Example 2 of Tillé and Rivest (N = 8 and 16), and the samples of the q x q grid balanced on rows, columns and anti-diagonals (q = 5, 7, 11), with the weights of the uniform design over them. |
| `time_preprocessing.py` | Section 6: time and memory to enumerate admissible directions and solve program (8), acceptance rates of rejection sampling, size of the stored minimax supports. |
| `make_tables.py` | Tables 1–5 and A1 from `results/`, including the exact shares of grid coset samples balanced on S2 and S4 and the exact main-effect weights of the coset designs. |
| `reproduce.sh` | Runs all of the above in order. |

## Data

| File | Content |
|---|---|
| `data/volcano61.csv` | Rows 14–74 (all 61 columns) of R's `datasets::volcano`; unit (i, j) has code 61 i + j. |
| `data/grid_y.csv`, `data/grid_X.csv` | Elevations in unit order, and the 0/1 indicators of rows, columns, anti-diagonals i+j and diagonals i−j (levels 1–60 of each; level 0 dropped). |
| `data/fact_p2_m7_Y.csv`, `data/fact_p3_m5_Y.csv` | The four responses of Section 7.1 (additive, mixed, interaction, noise), seeds 7 and 5. |
| `data/fact_p2_m7_X.csv`, `data/fact_p3_m5_X.csv` | Category indicators used as balancing variables by the cube method. |

## Results

Directions are stored as reduced row-echelon bases of the dual subspace D = V⊥ of the flat direction V.

| File | Content |
|---|---|
| `results/factorial_exact.csv` | Exact design effects of the coset designs of Section 7.1, with R² and adjusted R². |
| `results/factorial_minimax_mixtures.json` | Support (dual bases) and probabilities of the minimax mixtures of Section 7.1. |
| `results/factorial_lattice.csv` | Lattice sampling: strength, dual basis of the chosen array, exact design effects and worst case. |
| `results/factorial_designs.json` | Size of the support, worst case relative to the bound, positivity of second-order inclusion probabilities. |
| `results/factorial_cube_p2_m7.csv`, `results/factorial_cube_p3_m5.csv` | Monte Carlo summaries for the cube method (mean squared error, standard error, share exactly balanced, the main-effect, residual and cross parts of the mean squared error with their standard errors, and the main-effect weight). |
| `results/cube_samples/factorial_*.txt.gz`, `results/cube_samples/transcription_*.txt.gz` | Selected units of every factorial draw: package, and transcription with both landings. |
| `results/cube_factorial.csv` | Transcription on the factorial frames: design effects, shares exactly balanced, shares ending at an integral vertex, main-effect weights. |
| `results/replicated_populations.csv`, `results/replicated_populations.npz` | Quantiles of design effects and of ratios to the minimax design over 500 populations, and all design effects. |
| `results/attainment.json` | Enumeration of balanced samples where (8) stays above the bound; profiles and exact linear program for OA(16,7,2,2). |
| `results/varest_exact.csv` | Exact moments and coverage of the variance estimators. |
| `results/minimax_lp.csv` | Value of program (8) for every configuration solved (Table 1 shows a subset). |
| `results/lp_certificates.json` | Rational certificates: mixtures attaining the bound, and dual solutions with the exact minimax value where it is not attained. |
| `results/grid_exact.json` | Exact design effects of the grid designs, the evaluation of the variance estimator, and the sensitivity to ties among lattice directions. |
| `results/grid_lattice_ranking.csv` | The 58 directions admissible for S4 ranked by the length of the shortest dual vector; ties broken by increasing slope. |
| `results/grid_spectrum.json` | Residual Fourier energy by wrapped frequency length and on the dual lines of the chosen designs; seams of the wrapped grid. |
| `results/grid_balance_shares.json` | Exact share of the samples of each grid coset design that are balanced on S2 and on S4. |
| `results/mc_grid/*.csv` | Monte Carlo summaries for the cube, local pivotal and local cube methods on the grid. |
| `results/mc_grid_v2/*.csv` | Reruns with main-effect weights, and the local pivotal method with wrapped distances. |
| `results/cube_audit_S2.csv`, `results/cube_audit_S2_draws.csv` | Audit of the cube method balanced on rows and columns: summary, and one line per draw with the selected units. |
| `results/cube_audit_S2_decomposition.csv` | Decomposition of the mean squared error of the audited designs into main-effect, residual and cross parts. |
| `results/exhaustive_counts.json` | Exhaustive counts of Sections 1 and 6: vertices of the Example 2 polytopes, balanced samples of the q x q grids, flats among them and spectral weights of the uniform design. |
| `results/transcription_check.csv` | For each sample size and seed, whether the transcription returned the package's sample, and how many uniform numbers it used. |
| `results/preprocessing_costs.csv` | Preprocessing cost of the coset designs of Table 1 and Section 7.1 (one laptop core). |
| `results/table1.md` – `results/table5.md`, `results/tableA1.md` | The tables of the manuscript. |
