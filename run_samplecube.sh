#!/usr/bin/env bash
# samplecube of the package sampling (2.11) as a second implementation of the cube method (review round 11):
# factorial frames with both landings, 20 000 draws per configuration, and, for the cube method of BalancedSampling
# and both landings of samplecube, with the balancing variables in reverse order; grid on S2 (200 draws per sample size) and
# S4 (100 draws per sample size), in chunks with their own seeds, stored samples. samplecube is written in R and
# takes about 20 s per draw on S2 and 190 s on S4 on one core. Jobs whose summary file exists are skipped, so an
# interrupted run can be resumed.
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p results/mc_grid_scube results/cube_samples logs
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
{
for impl in scube1 scube2; do
  s=$([ $impl = scube1 ] && echo 21 || echo 23)
  [ -f results/factorial_${impl}_p2_m7.csv ] || echo "Rscript mc_cube.R data/fact_p2_m7_X.csv data/fact_p2_m7_Y.csv 2 7 8,16,32 20000 results/factorial_${impl}_p2_m7.csv $s results/cube_samples/factorial_${impl}_p2_m7.txt.gz $impl > logs/${impl}_p2_m7.log 2>&1"
  [ -f results/factorial_${impl}_p3_m5.csv ] || echo "Rscript mc_cube.R data/fact_p3_m5_X.csv data/fact_p3_m5_Y.csv 3 5 9,27,81 20000 results/factorial_${impl}_p3_m5.csv $((s+1)) results/cube_samples/factorial_${impl}_p3_m5.txt.gz $impl > logs/${impl}_p3_m5.log 2>&1"
done
for impl in bs scube1 scube2; do                        # the same with the balancing variables in reverse order
  s=$([ $impl = bs ] && echo 31 || ([ $impl = scube1 ] && echo 33 || echo 35))
  [ -f results/factorial_${impl}_rev_p2_m7.csv ] || echo "Rscript mc_cube.R data/fact_p2_m7_X.csv data/fact_p2_m7_Y.csv 2 7 8,16,32 20000 results/factorial_${impl}_rev_p2_m7.csv $s results/cube_samples/factorial_${impl}_rev_p2_m7.txt.gz $impl reversed > logs/${impl}_rev_p2_m7.log 2>&1"
  [ -f results/factorial_${impl}_rev_p3_m5.csv ] || echo "Rscript mc_cube.R data/fact_p3_m5_X.csv data/fact_p3_m5_Y.csv 3 5 9,27,81 20000 results/factorial_${impl}_rev_p3_m5.csv $((s+1)) results/cube_samples/factorial_${impl}_rev_p3_m5.txt.gz $impl reversed > logs/${impl}_rev_p3_m5.log 2>&1"
done
for k in 61 183; do
  for c in 0 1 2 3; do
    [ -f results/mc_grid_scube/scubeS4_k${k}_c$c.csv ] || echo "Rscript mc_grid.R $k 25 results/mc_grid_scube/scubeS4_k${k}_c$c.csv $((7000+k+c)) scubeS4 results/cube_samples/grid_scubeS4_k${k}_c$c.txt.gz > logs/scubeS4_k${k}_c$c.log 2>&1"
  done
  for c in 0 1; do
    [ -f results/mc_grid_scube/scubeS2_k${k}_c$c.csv ] || echo "Rscript mc_grid.R $k 100 results/mc_grid_scube/scubeS2_k${k}_c$c.csv $((6000+k+c)) scubeS2 results/cube_samples/grid_scubeS2_k${k}_c$c.txt.gz > logs/scubeS2_k${k}_c$c.log 2>&1"
  done
done
} | xargs -P "${JOBS:-10}" -I{} bash -c "{}"
