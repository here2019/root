#!/usr/bin/env bash
set -euo pipefail

echo "CPU:"
nproc || true
lscpu | head -n 30 || true

echo "Build:"
g++ -O3 -march=native -std=c++20 wheel105_pi.cpp -o wheel105_pi

echo "Run:"
for N in 100000000 200000000 500000000; do
  echo "=== N=$N ==="
  /usr/bin/time -f "time=%e s" ./wheel105_pi $N 16
done
