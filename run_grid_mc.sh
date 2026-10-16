#!/usr/bin/env bash
# Grid Monte Carlo runs used in the manuscript (seeds and replicate numbers as reported).
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p results/mc_grid logs
for k in 61 183; do
  Rscript mc_grid.R $k 400   results/mc_grid/cubeS4_k$k.csv  $((100+k)) cubeS4  > logs/cubeS4_k$k.log  2>&1 &
  Rscript mc_grid.R $k 500   results/mc_grid/lcubeS4_k$k.csv $((200+k)) lcubeS4 > logs/lcubeS4_k$k.log 2>&1 &
  Rscript mc_grid.R $k 2000  results/mc_grid/cubeS2_k$k.csv  $((300+k)) cubeS2  > logs/cubeS2_k$k.log  2>&1 &
  Rscript mc_grid.R $k 20000 results/mc_grid/lpm2_k$k.csv    $((400+k)) lpm2    > logs/lpm2_k$k.log    2>&1 &
done
wait
