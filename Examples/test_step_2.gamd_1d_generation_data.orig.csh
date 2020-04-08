#!/bin/csh

../2_gamd_reweight_run1d.py \
  -inp  b_orig_data/step_2.1d_gamd_generation_data.orig.list \
  -col  2       \
  -gamd fgf21-v16-s2.all_gamd.weights.dat.bz2 \
  -Emax 5       \
  -pwd  "`pwd`/b_orig_data" \
  -bin  0.5     \
  -temp 310     
