#!/usr/bin/env bash
# Stored Monte Carlo samples for the multi-surface test (seeds fixed in advance).
set -euo pipefail
cd "$(dirname "$0")/.."
for k in 61 183; do
  Rscript scripts/draw_mc.R lcubeS2 $k 2000 $((1000+k)) results/samples_lcubeS2_k$k.txt.gz > logs/lcubeS2_k$k.log 2>&1 &
  Rscript scripts/draw_mc.R cubeS2  $k 1000 $((2000+k)) results/samples_cubeS2_k$k.txt.gz  > logs/cubeS2_k$k.log  2>&1 &
done
wait
for k in 61 183; do
  Rscript scripts/draw_mc.R lpm2 $k 20000 $((3000+k)) results/samples_lpm2_k$k.txt.gz > logs/lpm2_k$k.log 2>&1
done
