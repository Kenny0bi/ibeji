#!/usr/bin/env bash
# Train every ibeji model set, one after another, resumable.
# Each set uses 4 workers; chromosomes with existing output are skipped.
set -euo pipefail
cd "$(dirname "$0")/.."

# Order: the cross-ancestry pair (YRI87, EUR87_r1), then EUR87_r2 for the within-ancestry
# noise floor, then EUR358 for the sample-size comparison, then the remaining draws.
SETS=(YRI87 EUR87_r1 EUR87_r2 EUR358 EUR87_r3 EUR87_r4 EUR87_r5)
CORES="${CORES:-4}"   # one worker per physical core; hyperthreads add little for glmnet

for set in "${SETS[@]}"; do
  echo "=== $(date '+%H:%M:%S') start $set"
  Rscript src/03_train_grex.R "$set" 1-22 "$CORES" 2>&1 | grep -v "built under R version" | tee -a "logs/train_${set}.log"
  echo "=== $(date '+%H:%M:%S') done $set"
done
