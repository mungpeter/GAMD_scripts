#!/bin/csh

../2_gamd_reweight_run2d.py \
  -inp  b_data/step_2.2d_gamd_generation_data.list \
  -pwd  "`pwd`" \
  -dir  b_data  \
  -col  2       \
  -gamd fgf21-v16-s2.all_gamd.weights.half.dat.bz2 \
  -Emax 5       \
  -job  amdweight_CE \
  -bin  0.5     \
  -Pmax 5       \
  -temp 310     
