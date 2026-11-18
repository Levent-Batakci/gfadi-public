#!/bin/bash
# The runs behind Fig. 1 (model problem of Section 1.1): one MATLAB process per (N, method), one at a
# time, appending to teaser/teaser.csv.  g-fADI runs first at each N, since the Krylov methods stop at
# the residual it reached.  Needs bash (Linux, macOS, or Git Bash on Windows) with matlab on the PATH.
# Usage: teaser/drive.sh 16 24 32 40 48     (rows already in the CSV are not skipped)
# Then draw the figure with make_teaser_fig in MATLAB.
cd "$(dirname "$0")"
out=teaser.csv
[ -f $out ] || echo "N,m,method,time,peakGB,its,res,projected,note" > $out
for N in "$@"; do
  for meth in gfadi fadi eksm arnoldi bs; do
    timeout 3600 matlab -batch "teaser_run($N, '$meth', '$out')" 2>&1 | grep -E "^N=|Error|error" 
  done
done
