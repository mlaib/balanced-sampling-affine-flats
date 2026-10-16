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
if [ "${QUICK:-0}" != "1" ]; then
  Rscript mc_cube.R data/fact_p2_m7_X.csv data/fact_p2_m7_Y.csv 2 7 8,16,32 20000 results/factorial_cube_p2_m7.csv 11
  Rscript mc_cube.R data/fact_p3_m5_X.csv data/fact_p3_m5_Y.csv 3 5 9,27,81 20000 results/factorial_cube_p3_m5.csv 12
  bash run_grid_mc.sh                       # Section 7.2: cube, local pivotal and local cube methods
  python3 audit_cube.py 2000 12             # Section 7.2: audit of the cube method on rows and columns
fi
python3 make_tables.py                      # Tables 1-4 (results/table1.md ... table4.md)
