#!/usr/bin/env bash
# Grid Monte Carlo reruns for review round 8: every design of Table 3 with all its draws and the original seeds
# (so the same draws as run_grid_mc.sh), recording the sample size, the Horvitz-Thompson estimate with the
# fixed denominator k and the main-effect weights of every draw; the cube and local cube samples are stored.
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p results/mc_grid_v3 results/cube_samples logs
for k in 61 183; do
  Rscript mc_grid.R $k 400   results/mc_grid_v3/cubeS4_k$k.csv    $((100+k)) cubeS4    results/cube_samples/grid_cubeS4_k$k.txt.gz  > logs/v3_cubeS4_k$k.log  2>&1 &
  Rscript mc_grid.R $k 500   results/mc_grid_v3/lcubeS4_k$k.csv   $((200+k)) lcubeS4   results/cube_samples/grid_lcubeS4_k$k.txt.gz > logs/v3_lcubeS4_k$k.log 2>&1 &
  Rscript mc_grid.R $k 2000  results/mc_grid_v3/cubeS2_k$k.csv    $((300+k)) cubeS2    results/cube_samples/grid_cubeS2_k$k.txt.gz  > logs/v3_cubeS2_k$k.log  2>&1 &
  Rscript mc_grid.R $k 20000 results/mc_grid_v3/lpm2_k$k.csv      $((400+k)) lpm2      > logs/v3_lpm2_k$k.log      2>&1 &
  Rscript mc_grid.R $k 20000 results/mc_grid_v3/lpm2torus_k$k.csv $((500+k)) lpm2torus > logs/v3_lpm2torus_k$k.log 2>&1 &
done
wait
