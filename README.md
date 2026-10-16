# Exactly balanced sampling on overlapping stratifications with affine flats

Code, data and results for the manuscript of the same title by Mohamed Laib. Every number in the
tables of the manuscript is produced by the scripts below from the files in `data/`, and every
intermediate result is stored in `results/`.

## Requirements

- Python 3.12 with numpy 2.4.6 and scipy 1.18.0 (`requirements.txt`), for all exact results.
- R 4.3.3 with BalancedSampling 2.1.1, for the Monte Carlo runs of the cube, local pivotal and local cube methods.

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
| `test_theory.py` | Brute-force checks (support enumeration) of every exact statement in Sections 2–6; prints the number of checks per statement. |
| `study_factorial.py` | Section 7.1 (complete factorial frames 2^7 and 3^5): responses, exact design effects, minimax mixtures; writes the inputs of `mc_cube.R`. |
| `study_lattice.py` | Section 7.1: lattice sampling with a linear orthogonal array of maximal strength (Table 2, lattice column). |
| `mc_cube.R` | Section 7.1: Monte Carlo for the cube method (random unit order before each draw). |
| `study_minimax_lp.py` | Table 1: minimax mixing distributions by linear programming, reduced by the monomial group. |
| `certify_lp.py` | Exact rational certificates for program (8): primal mixtures where the bound is attained, dual solutions where it is not. |
| `study_grid.py` | Section 7.2 (61 x 61 elevation grid): exact design effects, ranking of the lattice directions, exact evaluation of the variance estimator (Table 4), residual Fourier energy; writes the inputs of `mc_grid.R`. |
| `mc_grid.R` | Section 7.2: Monte Carlo for the cube, local pivotal (lpm2) and local cube (lcube) methods. |
| `run_grid_mc.sh` | Runs `mc_grid.R` for every design and sample size of Section 7.2 with the seeds and replicate numbers of Appendix A. |
| `audit_cube.py` | Section 7.2: Python transcription of the flight and landing phases of `cube` in BalancedSampling 2.1.1, instrumented, with the flight optionally continued until no move preserves all balancing equations; `--decompose` splits each mean squared error into main-effect, residual and cross parts. |
| `time_preprocessing.py` | Section 6: time and memory to enumerate admissible directions and solve program (8), acceptance rates of rejection sampling, size of the stored minimax supports. |
| `make_tables.py` | Tables 1–4 from `results/`, including the exact shares of grid coset samples balanced on S2 and S4. |
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
| `results/factorial_cube_p2_m7.csv`, `results/factorial_cube_p3_m5.csv` | Monte Carlo summaries for the cube method (mean squared error, standard error, share exactly balanced). |
| `results/minimax_lp.csv` | Value of program (8) for every configuration solved (Table 1 shows a subset). |
| `results/lp_certificates.json` | Rational certificates: mixtures attaining the bound, and dual solutions with the exact minimax value where it is not attained. |
| `results/grid_exact.json` | Exact design effects of the grid designs, the evaluation of the variance estimator, and the sensitivity to ties among lattice directions. |
| `results/grid_lattice_ranking.csv` | The 58 directions admissible for S4 ranked by the length of the shortest dual vector; ties broken by increasing slope. |
| `results/grid_spectrum.json` | Residual Fourier energy by wrapped frequency length and on the dual lines of the chosen designs; seams of the wrapped grid. |
| `results/grid_balance_shares.json` | Exact share of the samples of each grid coset design that are balanced on S2 and on S4. |
| `results/mc_grid/*.csv` | Monte Carlo summaries for the cube, local pivotal and local cube methods on the grid. |
| `results/cube_audit_S2.csv`, `results/cube_audit_S2_draws.csv` | Audit of the cube method balanced on rows and columns: summary, and one line per draw with the selected units. |
| `results/cube_audit_S2_decomposition.csv` | Decomposition of the mean squared error of the audited designs into main-effect, residual and cross parts. |
| `results/preprocessing_costs.csv` | Preprocessing cost of the coset designs of Table 1 and Section 7.1 (one laptop core). |
| `results/table1.md` – `results/table4.md` | The tables of the manuscript. |
