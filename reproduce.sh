#!/usr/bin/env bash
# Rebuild every result and table of the manuscript.  Exact computations take minutes; the Monte Carlo
# runs take hours (set QUICK=1 to skip them and rebuild the tables from the stored results).
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p logs
python3 test_theory.py                      # brute-force checks of Sections 2-4
python3 study_factorial.py                  # Section 5.1: exact design effects; writes data/fact_*
python3 study_lattice.py                    # Section 5.1: lattice sampling
python3 study_minimax_lp.py                 # Table 1
python3 certify_lp.py                       # exact certificates for Table 1 (attainment and non-attainment)
python3 study_grid.py                       # Section 5.2: exact results, lattice ranking, spectrum, inputs of Table 3
python3 time_preprocessing.py               # Section 4: preprocessing cost of the minimax mixtures
python3 count_balanced_samples.py           # Sections 4 and 5.1: exhaustive counts of balanced samples and vertices
python3 settle_attainment.py                # Section 4: no exactly balanced design attains (4) where (8) stays above it (needs gcc)
python3 crosscheck_oa16.py                  # Appendix A: independent check of the OA(16,7,2,2) enumeration (needs gcc)
python3 check_certificates.py               # Appendix A: standalone check of the certificates (standard library only)
python3 varest_exact.py                     # Section 4 and Supplement S6: exact evaluation of the variance estimators
python3 check_ranks.py                      # Appendix A: ranks of the regression estimator, by singular values and exactly
if [ "${QUICK:-0}" != "1" ]; then
  mkdir -p results/cube_samples
  Rscript mc_cube.R data/fact_p2_m7_X.csv data/fact_p2_m7_Y.csv 2 7 8,16,32 20000 results/factorial_cube_p2_m7.csv 11 results/cube_samples/factorial_p2_m7.txt.gz
  Rscript mc_cube.R data/fact_p3_m5_X.csv data/fact_p3_m5_Y.csv 3 5 9,27,81 20000 results/factorial_cube_p3_m5.csv 12 results/cube_samples/factorial_p3_m5.txt.gz
  python3 cube_factorial.py 20000 5         # Section 5.1: cube method with the flight continued to a vertex
  bash run_grid_mc_v3.sh                    # Section 5.2, Tables 3 and S5: every grid design with the original seeds (run_grid_mc.sh and
                                            # run_grid_mc_v2.sh, the earlier runs with the same seeds, are kept for the record)
  python3 audit_cube.py 2000 12             # Section 5.2 and Supplement S4: audit of the cube method on rows and columns
  python3 cube_grid_fixed.py 400 8          # Section 5.2: cube method on S4 with directions by QR factorisation
  python3 check_rref_defect.py 4            # Supplement S4: the size exceptions come from the package's row reduction
  Rscript check_size_drift.R                # Supplement S4: sizes returned by the package's cube and lcube on a 19 x 19 grid
  python3 validate_transcription.py         # Supplement S4: the transcription returns the package's samples given R's random numbers
  bash run_samplecube.sh                    # Section 5 and Supplement S2: samplecube of the package sampling, and the runs in reverse order
fi
python3 ht_factorial.py                     # Tables 2, S3, S4 and S6: Horvitz-Thompson mean of the cube variants from the stored samples
python3 replicate_populations.py 500        # Supplement S3, Table S1: 500 further populations (uses the stored cube samples)
python3 bootstrap_bank.py 200 500           # Supplement S3: Monte Carlo standard errors of the cube quantiles (resampled draws)
python3 samplecube_grid.py                  # Section 5.2: summaries of the grid runs of samplecube
python3 bound_centred.py                    # Section 5, Tables 3, S2 and S5: lower bounds (S1) with the size correction
python3 make_tables.py                      # Tables 1-3, A.1 and S1-S7 (results/table1.md ... table3.md, tableA.1.md, tableS1.md ... tableS7.md)
