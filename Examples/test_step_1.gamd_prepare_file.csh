#!/bin/csh

## this scirpt is only as an example and 
## won't be run in this /example folder
../1_gamd_reweight_prep.py \
   -mf 1_v16 2_v16a 3_v16r 4_v16ar \
   -sf 1_run 2_run                 \
   -log gamd.1.log.bz2 gamd.2.log.bz2 gamd.3.log.bz2 \
   -o fgf21_variants_gamd          \
   -t 310 -s 1

