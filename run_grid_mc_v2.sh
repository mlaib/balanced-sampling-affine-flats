#!/usr/bin/env bash
# Grid Monte Carlo reruns for review round 7: the same seeds as run_grid_mc.sh (so the same draws), with the
# main-effect weights t_M recorded; the slow cube (S4) and local cube runs use their first 100 draws only.
# Adds the local pivotal method with wrapped distances (lpm2torus, seeds 500 + k).
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p results/mc_grid_v2 logs
for k in 61 183; do
  Rscript mc_grid.R $k 100   results/mc_grid_v2/cubeS4_k$k.csv    $((100+k)) cubeS4    > logs/v2_cubeS4_k$k.log    2>&1 &
  Rscript mc_grid.R $k 100   results/mc_grid_v2/lcubeS4_k$k.csv   $((200+k)) lcubeS4   > logs/v2_lcubeS4_k$k.log   2>&1 &
  Rscript mc_grid.R $k 2000  results/mc_grid_v2/cubeS2_k$k.csv    $((300+k)) cubeS2    > logs/v2_cubeS2_k$k.log    2>&1 &
  Rscript mc_grid.R $k 20000 results/mc_grid_v2/lpm2_k$k.csv      $((400+k)) lpm2      > logs/v2_lpm2_k$k.log      2>&1 &
  Rscript mc_grid.R $k 20000 results/mc_grid_v2/lpm2torus_k$k.csv $((500+k)) lpm2torus > logs/v2_lpm2torus_k$k.log 2>&1 &
done
wait
