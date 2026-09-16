#!/usr/bin/env bash
# The remaining within-ancestry noise-floor pairs, alongside src/auto_downstream.sh.
#
# The noise floor is the control the whole decomposition rests on: it is what shows
# that most of the apparent cross-ancestry weights difference is training noise. That
# claim currently stands on a single pair of European draws. With five draws there are
# ten possible pairs and each decomposition takes about a minute, so the floor can be
# a distribution instead of one number.
#
# auto_downstream.sh already runs r1-r2, r1-r3 and r2-r3. This runs the other seven.
# The step names are disjoint from that script's, so the two loops never race for the
# same marker, and each pair waits until both of its training sets have all 22
# chromosomes. Markers in results/.done/ mean a restart never repeats finished work.
cd "$(dirname "$0")/.."
DONE=results/.done
mkdir -p "$DONE" logs
LOG=logs/auto_extra_pairs.log
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
  if nice -n 10 "$@" > "logs/auto_${name}.log" 2>&1; then
    touch "$DONE/$name"
    echo "$(date '+%a %H:%M') done  $name" >> "$LOG"
  else
    echo $((tries + 1)) > "$tries_file"
    echo "$(date '+%a %H:%M') FAILED $name (attempt $((tries + 1)) of 3, see logs/auto_${name}.log)" >> "$LOG"
  fi
}

PAIRS=("EUR87_r1 EUR87_r4" "EUR87_r1 EUR87_r5" "EUR87_r2 EUR87_r4" "EUR87_r2 EUR87_r5" \
       "EUR87_r3 EUR87_r4" "EUR87_r3 EUR87_r5" "EUR87_r4 EUR87_r5")
echo "$(date '+%a %H:%M') auto_extra_pairs started, ${#PAIRS[@]} pairs queued" >> "$LOG"

while true; do
  remaining=0
  for pair in "${PAIRS[@]}"; do
    set -- $pair
    name="decomp_${1}_vs_${2}"
    if [ ! -f "$DONE/$name" ]; then
      remaining=$((remaining + 1))
      complete "$1" && complete "$2" && step "$name" Rscript src/05_decompose.R "$1" "$2" "$WORKERS"
    fi
  done
  if [ "$remaining" -eq 0 ]; then
    echo "$(date '+%a %H:%M') all extra pairs complete; exiting" >> "$LOG"
    break
  fi
  sleep 300
done
