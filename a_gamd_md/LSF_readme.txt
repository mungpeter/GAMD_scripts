
## function of each script here

XXX = peptide (<15aa), kinase (~200-300aa), klotho (~1000aa)

# /explicit_input_XXX
  the folder contains input files for AMBER MD simulation steps, from minimization, 
  heating, equilibration, to GAMD initial parameter collection to GAMD production 
  run

# LSF_amber_gamd.XXX-min
  Batch scheduler submitting jobs to LSF quening system. It is written like a shell 
  script. This one runs specifically the initial system minimization/equilibration/
  gamd_initiation.

# rerun_amber.csh
  At end of LSF...-min script running, it will call this script to loop over the
  gamd production run script LSF_...-run, e.g. looping over the production run of
  50ns 5x times will end up running the production script 5x and get 250ns final
  trajectory. The number of looping is pre-defined in LSF...-min.

# LSF_amber_gamd.XXX-run
  LSF script to do gamd_production run, which can be looped over.

## All these items should be in the same directory level unless changed in the LSF
   scripts.
