#!/usr/bin/env bash
# Rebuild every result and table of the manuscript.  Exact computations take minutes; the Monte Carlo
# runs take hours (set QUICK=1 to skip them and rebuild the tables from the stored results).
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p logs
python3 test_theory.py                      # brute-force checks of Sections 2-6
python3 study_factorial.py                  # Section 7.1: exact design effects; writes data/fact_*
python3 study_lattice.py                    # Section 7.1: lattice sampling
python3 study_minimax_lp.py                 # Table 1
python3 certify_lp.py                       # exact certificates for Table 1 (attainment and non-attainment)
python3 study_grid.py                       # Section 7.2: exact results, lattice ranking, spectrum, Table 4 inputs
python3 time_preprocessing.py               # Section 6: preprocessing cost of the minimax mixtures
python3 count_balanced_samples.py           # Sections 6 and 7.1: exhaustive counts of balanced samples and vertices
python3 settle_attainment.py                # Section 6: no exactly balanced design attains (4) where (8) stays above it (needs gcc)
python3 varest_exact.py                     # Section 7.3: exact evaluation of the variance estimators
if [ "${QUICK:-0}" != "1" ]; then
  mkdir -p results/cube_samples
  Rscript mc_cube.R data/fact_p2_m7_X.csv data/fact_p2_m7_Y.csv 2 7 8,16,32 20000 results/factorial_cube_p2_m7.csv 11 results/cube_samples/factorial_p2_m7.txt.gz
  Rscript mc_cube.R data/fact_p3_m5_X.csv data/fact_p3_m5_Y.csv 3 5 9,27,81 20000 results/factorial_cube_p3_m5.csv 12 results/cube_samples/factorial_p3_m5.txt.gz
  python3 cube_factorial.py 20000 5         # Section 7.1: cube method with the flight continued to a vertex
  bash run_grid_mc.sh                       # Section 7.2: cube, local pivotal and local cube methods
  bash run_grid_mc_v2.sh                    # Section 7.2: main-effect weights, local pivotal method with wrapped distances
  python3 audit_cube.py 2000 12             # Section 7.2: audit of the cube method on rows and columns (and its MSE decomposition)
  python3 validate_transcription.py         # Appendix A: the transcription returns the package's samples given R's random numbers
fi
python3 replicate_populations.py 500        # Appendix A, Table A1: 500 further populations (uses the stored cube samples)
python3 make_tables.py                      # Tables 1-5 and A1 (results/table1.md ... table5.md, tableA1.md)
