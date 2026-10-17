#!/usr/bin/env bash
# Rebuild every result and table of the manuscript.  Exact computations take minutes; the Monte Carlo
# runs take hours (set QUICK=1 to skip them and rebuild the tables from the stored results).
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p logs
python3 test_theory.py                      # brute-force checks of Sections 2-6
python3 study_factorial.py                  # Section 7.1: exact design effects; writes data/fact_*
python3 study_lattice.py                    # Section 7.1: lattice sampling
python3 study_minimax_lp.py                 # Table 2
python3 certify_lp.py                       # exact certificates for Table 2 (attainment and non-attainment)
python3 study_grid.py                       # Section 7.2: exact results, lattice ranking, spectrum, inputs of Tables 5 and 6
python3 time_preprocessing.py               # Section 6: preprocessing cost of the minimax mixtures
python3 count_balanced_samples.py           # Sections 6 and 7.1: exhaustive counts of balanced samples and vertices
python3 settle_attainment.py                # Section 6: no exactly balanced design attains (4) where (8) stays above it (needs gcc)
python3 crosscheck_oa16.py                  # Appendix A: independent check of the OA(16,7,2,2) enumeration (needs gcc)
python3 check_certificates.py               # Appendix A: standalone check of the certificates (standard library only)
python3 varest_exact.py                     # Section 7.3: exact evaluation of the variance estimators
if [ "${QUICK:-0}" != "1" ]; then
  mkdir -p results/cube_samples
  Rscript mc_cube.R data/fact_p2_m7_X.csv data/fact_p2_m7_Y.csv 2 7 8,16,32 20000 results/factorial_cube_p2_m7.csv 11 results/cube_samples/factorial_p2_m7.txt.gz
  Rscript mc_cube.R data/fact_p3_m5_X.csv data/fact_p3_m5_Y.csv 3 5 9,27,81 20000 results/factorial_cube_p3_m5.csv 12 results/cube_samples/factorial_p3_m5.txt.gz
  python3 cube_factorial.py 20000 5         # Section 7.1: cube method with the flight continued to a vertex
  bash run_grid_mc_v3.sh                    # Section 7.2, Tables 5 and 6: every grid design with the original seeds (run_grid_mc.sh and
                                            # run_grid_mc_v2.sh, the earlier runs with the same seeds, are kept for the record)
  python3 audit_cube.py 2000 12             # Section 7.2: audit of the cube method on rows and columns (and its MSE decomposition)
  python3 cube_grid_fixed.py 400 8          # Section 7.2: cube method on S4 with directions by QR factorisation
  python3 check_rref_defect.py 4            # Supplement S2: the size exceptions come from the package's row reduction
  Rscript check_size_drift.R                # Supplement S2: sizes returned by the package's cube and lcube on a 19 x 19 grid
  python3 validate_transcription.py         # Supplement S2: the transcription returns the package's samples given R's random numbers
fi
python3 ht_factorial.py                     # Tables 3 and 4: Horvitz-Thompson mean of the cube variants from the stored samples
python3 replicate_populations.py 500        # Supplement S3, Table S1: 500 further populations (uses the stored cube samples)
python3 make_tables.py                      # Tables 2-7 and S1 (results/table2.md ... table7.md, tableS1.md)
