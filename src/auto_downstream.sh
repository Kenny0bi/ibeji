#!/usr/bin/env bash
# Run every downstream analysis as soon as its inputs exist, with no manual step.
#
# Runs beside src/run_all_training.sh. Every 5 minutes it checks which training sets are
# complete (22 chromosome weight files) and runs, once each:
#   - TWAS for every complete set
#   - decompositions for every pair whose two sets are complete
#   - model QC once the core sets are done
# Finished steps leave a marker in results/.done/, so a restart never repeats work.
# Failed steps are logged as FAILED and retried on the next pass (up to 3 attempts).
cd "$(dirname "$0")/.."
DONE=results/.done
mkdir -p "$DONE" logs
LOG=logs/auto_downstream.log
WORKERS=2

complete() {
  [ "$(ls results/models/"$1"/chr*.weights.tsv 2>/dev/null | wc -l | tr -d ' ')" -eq 22 ]
}

step() {
  local name=$1; shift
  [ -f "$DONE/$name" ] && return 0
  local tries_file="$DONE/$name.tries"
  local tries=$(cat "$tries_file" 2>/dev/null || echo 0)
  [ "$tries" -ge 3 ] && return 0
  echo "$(date '+%a %H:%M') start $name" >> "$LOG"
  if "$@" > "logs/auto_${name}.log" 2>&1; then
    touch "$DONE/$name"
    echo "$(date '+%a %H:%M') done  $name" >> "$LOG"
  else
    echo $((tries + 1)) > "$tries_file"
    echo "$(date '+%a %H:%M') FAILED $name (attempt $((tries + 1)) of 3, see logs/auto_${name}.log)" >> "$LOG"
  fi
}

SETS=(YRI87 EUR87_r1 EUR87_r2 EUR358 EUR87_r3 EUR87_r4 EUR87_r5)
echo "$(date '+%a %H:%M') auto_downstream started" >> "$LOG"

while true; do
  for s in "${SETS[@]}"; do
    complete "$s" && step "twas_$s" Rscript src/06_twas.R "$s" "$WORKERS"
  done

  if complete YRI87; then
    for e in EUR87_r1 EUR87_r2 EUR358 EUR87_r3 EUR87_r4 EUR87_r5; do
      complete "$e" && step "decomp_${e}_vs_YRI87" Rscript src/05_decompose.R "$e" YRI87 "$WORKERS"
    done
  fi
  for pair in "EUR87_r1 EUR87_r2" "EUR87_r1 EUR87_r3" "EUR87_r2 EUR87_r3"; do
    set -- $pair
    complete "$1" && complete "$2" && step "decomp_${1}_vs_${2}" Rscript src/05_decompose.R "$1" "$2" "$WORKERS"
  done

  if complete YRI87 && complete EUR87_r1 && complete EUR87_r2 && complete EUR358; then
    step model_qc_core Rscript src/08_model_qc.R
  fi

  all_done=1
  for s in "${SETS[@]}"; do complete "$s" || all_done=0; done
  if [ "$all_done" -eq 1 ]; then
    step model_qc_all Rscript src/08_model_qc.R
    echo "$(date '+%a %H:%M') all training sets complete; downstream loop exiting" >> "$LOG"
    break
  fi
  sleep 300
done
