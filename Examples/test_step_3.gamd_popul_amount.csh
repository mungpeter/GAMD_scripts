#!/bin/csh

## the script will use glob.glob to find the full path of the files
## use -pwd if path-name has characters non-compatible with Linux (' ', '\', '(', etc)

../3_gamd_extract_popul.py \
  -in   b_orig_data/fgf21-v16-s2-cp.all.204-208.rmsd.txt.bz2 b_orig_data/fgf21-v16-s2-cp.all.dist_p205-w208.txt.bz2 \
  -col  2     \
  -pwd  "`pwd`" \
  -list b_orig_data/step_3.gamd_popul_amount.list
