## GAMD_scripts
**Scripts to run _Gaussian Accelerated MD_ (**GAMD**) simulations and post-processing**

```
  author: Peter M.U. Ung @ MSSM/Yale
  vers:   1.    2019.05
```
See the original webpage and papers for more details:

[Hyperdyanmics: Unconstrained Gaussian Accelerated MD simulation method](https://gamd.ucsd.edu/)

_Reference 1_: [Miao, Y.; Feher, V. A.; McCammon, J. A., Gaussian Accelerated Molecular Dynamics: Unconstrained Enhanced Sampling and Free Energy Calculation. J. Chem. Theory Comput. 2015, 11, 3584-3595.](https://doi.org/10.1021/acs.jctc.5b00436)

_Reference 2_: [Y. Miao and J. A. McCammon (2016), Graded activation and free energy landscapes of a muscarinic G protein-coupled receptor. Proc Natl Acad Sci U S A, 113 (43): 12162–12167.](https://doi.org/10.1073/pnas.1614538113)

######################################################################################
- **Example scripts and AMBER MD setup to run GAMD**
```
/a_gamd_md
  - LSF_amber_gamd.XXX-min   [ LSF script to run AMBER GAMD simulation in the initial state, plus 4_gamd_init ]
  - LSF_amber_gamd.XXX-run   [ LSF script to run AMBER GAMD simulation in production state ]
  - rerun_amber.csh          [ rerun a production run from previous restart file with new number ]

- Required packages:
    AMBER       # stable: 16+, nvidia V100 CUDA-enabled
    cpptraj     # github version is better
```
- Contains AMBER scripts to run GAMD simulations. See the LSF scripts within for submission to HPC.


- Example data
```
/b_data
/b_orig_data
```
- Processed Data files of _GAMD Energy_ and _Structural Metrics_. The Energy and Metric files should have **matching number of data points**. 
*/b_orig_data* contains the _original_ dataset, **32,000** data points from 2 trajectories of 320ns GAMD run, taken at _20-ps interval_.
*/b_data* is the strip-down version of */b_orig_data*, **16,000** data points each, taken at _40-ps interval_.

- Example results
```
/c_result
/c_orig_result
```
- _2D heatmaps of free energy and conformational population landscapes image files_. Compare the two sets, you will see the strip-down version has less detail and some of the less-populated states are not shown. 
- This is also an example of the effect of **under-sampling** when generating the result. Recommend to take data from _no less than 20-ps interval and at least 500-ns_ worth of data for short peptide analysis. For larger system with large-scale movement, larger interval would be okay.

######################################################################################
- **Collect and Preparae the Weighted GAMD parameters from various GAMD trajectories**
```
> 1_gamd_reweight_prep.py
   -mf <a b c d ...>     [ parent folder(s) of GAMD data ]
   -sf <a b c d ...>     [ subfolder(s) of GAMD data ]
   -log <a b c d ...>    [ raw GAMD Energy log file(s) ]
   -o <prefix>           [ output prefix of processed GAMD Energy data ]
   
 Optional:
   -t <310>              [ simulation temperature K (def: 310) ]
   -s <1>                [ interval of data point to take (def: 1) ]
   -h, --help            [ show this help message and exit ]
   
 e.g.>  1_gamd_reweight_prep.py              \
            -mf  1_v16 2_v16a 3_v16r 4_v16ar \
            -sf  1_run 2_run                 \
            -log gamd.1.log.bz2 gamd.2.log.bz2 gamd.3.log.bz2 \
            -o  fgf21_vars_gamd_weighted     \
            -t  310                          \
            -s  1
```

- Preparation of the converted **GAMD Energy summary file** from individual raw GAMD Energy files found in _folders/subfolders_. Each raw file has header comments and the initial GAMD file (4_gamd_init) has a row of _0.00000000_ that needs to be removed. This script takes care of that. 
- This script also converts the GAMD energy from eV to kcal/mol at the Temperature used in the simulation, as instructed in the GAMD manual: _prepare_gamd_log.txt_


```
## For AMBER 14+
## shell script to convert data in gamd.log into formatted data
awk 'NR%1==0' gamd.log | awk '{print ($8+$7)/(0.0019872036*310)" " $2 " " ($8+$7)}' > gamd.weights.dat


- Running at this level of the directory structure*
  ----------1_v16____1_run_____gamd.1.log.bz2
      |      |         |_____gamd.2.log.bz2
      |      |         |_____gamd.3.log.bz2
      |      |
      |      |_____2_run_____gamd.1.log.bz2
      |                |_____ ...
      |
      |---2_v16a___1_run_____ ...
      |      |_____2_run_____ ...
      |
      |---3_v16r___ ...
      |---4_v16ar__ ...
```
- See the example file: **test_step_1.gamd_prepare_file.csh**

############
- **Analysis of GAMD trajectories to generate 1D energy results - estimation of energy well-depth**
```
> 2_gamd_reweight_run1d.py
  -inp  <input list>         [ list of input files ]
  -col  <column>             [ column number to be read ]
  -gamd <GAMD weight file>   [ GAMD weight file ]
  -Emax <Emax>               [ Max energy above 0 kcal/mol (def: 4) ]

Optional:
  -pwd  <path>               [ path to working directory, beware of characters' \()' (def: ./) ]
  -dir  <path>               [ path to data file, on top of -pwd (def: ./) ]
  
  -bin     <bin size>        [ General bin size for X-axes (def: 0.5) ]
  -Xdim    <Xmin Xmax>       [ Data-range in X-dimension (def: auto) ]
  -temp    <temperature>     [ Temperature K (def: 310) ]
  -histcut <histo cutoff>    [ Data cutoff (def: 25) ]

  -contour <contour>         [ Figure Contour per integer (def: 4) ]
  -smooth  <smooth>          [ Contour smoothening (def: None | Sug: 1.15) ]
  -dpi     <dpi>             [ Figure Resolution (def: 200) ]
  -img     <img>             [ Figure image format: png|svg|eps|pdf (def: png) ]

  -h, --help                 [ show this help message and exit ]

e.g.> ./2_gamd_reweight_run1d.py      \
  -inp  ./orig_data/step_2.1d_gamd_generation_data.orig.list \
  -gamd fgf21-v16-s2.all_gamd.weights.dat.bz2 \
  -col  2         \
  -Emax 5         \
                  \
  -pwd  /PATH-TO-WORKING-DIRECTORY \
  -dir  orig_data \
  -bin  0.5       \
  -temp 310
```

- This script generates the _GAMD result in **1D** figures_. It is only used to calculate and visual the _maximum depth of the energy well (**Emax**)_ for use in generating energy landscape in 2D analysis.
- See the example file: **test_step_2.gamd_1d_generation_data.orig.csh**


###########
- **Analysis of GAMD trajectories to generate 1D energy results - estimation of energy well-depth**
```
> 2_gamd_reweight_run2d.py
  -inp  <input list>         [ list of input files ]
  -col  <column>             [ column number to be read ]
  -gamd <GAMD weight file>   [ GAMD weight file ]
  -Emax <Emax>               [ Max energy above 0 kcal/mol (def: 4) ]

Optional:
  -pwd  <path>               [ path to working directory, beware of characters' \()' (def: ./) ]
  -dir  <path>               [ path to input data, on top of -pwd (def: ./) ]
  
  -job  <reweight method>    [ Reweighting method: noweight,  weighthist,   amd_time, amd_dV, 
                                                   amdweight, amdweight_MC, amdweight_CE
                               (Def: amdweight_CE) ]
  -temp <temperature>        [ Temperature K (def: 310) ]

  -bin   <bin size>          [ General bin size for both X/Y axes (def: 0.5) ]
  -Xbin  <Xbin size>         [ Bin size on X-axis, supersede -bin (def: auto) ]
  -Ybin  <bin size>          [ Bin size on Y-axis, supersede -bin (def: auto) ]
  -Xdim  <Xmin Xmax>         [ Data range of X-dimension (def: auto) ]
  -Ydim  <Ymin Ymax>         [ Data range of Y-dimension (def: auto) ]
  -Pmax  <Pmax>              [ Max Population percentage (def: 5) ]
  -order <MC order>          [ Order of MC formula (def: 10) ]
  -fit   <Fit switch>        [ Fitting data (default: False) ]

  -histcut <histo cutoff>    [ Data cutoff (def: 25) ]
  -contour <contour>         [ Figure Contour per integer step (def: 4) ]
  -smooth  <smooth>          [ Contour smoothening (def: None | Sug: 1.15) ]
  -dpi     <dpi>             [ Figure Resolution (def: 200) ]
  -img     <img>             [ Figure image format: png|svg|eps|pdf (def: png) ]

  -h, --help                 [ show this help message and exit ]
  
e.g.> ./2_gamd_reweight_run2d.py       \
  -inp  ./orig_data/step_2.2d_gamd_generation_data.orig.list \
  -gamd fgf21-v16-s2.all_gamd.weights.dat.bz2 \
  -col  2            \
  -Emax 4            \
                     \
  -pwd  /PATH-TO-WORKING-DIRECTORY \
  -dir  orig_data    \
  -job  amdweight_CE \
  -bin 0.5           \
  -Xbin 0.25         \
  -Pmax 5            \
  -temp 310          \
  -smooth 1.2        \
  -contour 4
```
- This script is the **MAIN** script to summarize _GAMD results into **2D** figures_ of **free energy** and **conformational population** landscapes. It takes the prepared GAMD Energy summary file and Structural Metrics to generate landscape files.
- See the example file: **test_step_2.gamd_2d_generation_data.orig.csh**

**2D energy landscape plot to identify potential energy wells and activation energy barrier sampled in GAMD runs.**

- One key thing about the conformational free energy landscape (measured in potential-of-mean force) is that GAMD has to have accessed that particular conformation in order to get an energy value. If GAMD couldn't sample a particular space, either due to undersampling from "too-short" simulations, or it is energetically inaccessable for GAMD, then there would be a gapping "high-energy" hole in the middle of a low-energy well. In that case, either extend the simulation, or change the size of '-bin' so that with averaging of a different size of group, the missing "hole" can be removed.

![2D energy landscape plot to identify potential energy wells and activation energy barrier sampled in GAMD runs](https://github.com/mungpeter/GAMD_scripts/blob/master/Examples/c_result/2D_dG_surf.fgf21-v16-s2-cp.all-half.204-208.rmsd.dist_s204-w207.dG.png)

**2D Conformational density plot to identify clusters of population actually sampled frequently in GAMD runs.**

- This is a direct representation of the GAMD sampling. In cases where the free-energy minima are accessible in the simulation conditions, the free-energy minima should correspond to the population density. This indicates that the measured region is conformationally well-behaved. However, if the free-energy minima is not really accessible due to either conformational constraints or lowe frequency of accessing due to high entropy, then the free-energy landscape and conformational density plots will be different. This is a tell-tale sign that the measured region has high entropy/flexibility.

![2D Population density plot to identify clusters of population actually sampled frequently in GAMD runs](https://github.com/mungpeter/GAMD_scripts/blob/master/Examples/c_result/2D_dG_surf.fgf21-v16-s2-cp.all-half.204-208.rmsd.dist_s204-w207.popul.png)

- When the x- and y-axis parameters were chosen that correctly reflect the energetically relevant conformational changes, the free-energy landscape (potential-of-mean force plot, above) and conformational population landscape plot (below) would see well defined clusters of population.

###########
- **Extract the quantitative number of population in the designated area in the density map**
```
> 3_gamd_extract_popul.py
  -in <X-input  Y-input>  [ Input data files ]
  -col <column>           [ Column number in file to read ]
  -list <range list>      [ List of X- and Y-ranges to collect, x-/y-ranges
                            separated by delimiters ',|;|:' ]
                            
Optional:
  -pwd <work directory>   [ Path to working directory, beware of characters' \()' (def: ./) ]
  -skip <skip line>       [ Skip every other X input lines (def: 0) ]
  -xrange <xmin xmax>     [ Range of X-dimension to collect ] * not use if -list
  -yrange <ymin ymax>     [ Range of Y-dimension to collect ] * not use if -list

  -h, --help              [ show this help message and exit ]

#####
# example input for -list file: 
  # 0 1 , 3 5.5         # population_1: x-range [0, 1] ;    y-range [3, 5.5]
  # 3.5 6.0 ; 3 5.5     # population_2: x-range [3.5,6.0] ; y-range [3, 5.5]
  # 3.5 5.0 : 10 15     # population_2: x-range [3.5,5.0] ; y-range [10, 15]
#####  

e.g.> ./3_gamd_extract_popul.py        \
  -in ./data/fgf21-v16-s2-cp.all.204-208.rmsd.txt.bz2 ./data/fgf21-v16-s2-cp.all.dist_p205-w208.txt.bz2 \
  -col 2                          \
  -pwd /PATH-TO-WORKING-DIRECTORY \
  -list ./data/step_3.gamd_popul_amount.list

-------------------------------
# X-range [ 0.0 , 1.0 ] | Y-range [ 0.0 , 4.0 ]
 % Population: 36.15

# X-range [ 1.0 , 3.0 ] | Y-range [ 0.0 , 4.0 ]
 % Population: 20.29

# X-range [ 1.0 , 4.0 ] | Y-range [ 4.0 , 10.0 ]
 % Population: 40.53

# % Total Population: 96.98

```

- This script counts the _frequency MD frames falling within specified ranges in a 2D landscape_. This is essentially binning to see what is the size of different conformational populations.
- See the example file: **test_step_3.gamd_popul_amount.csh**


Population number extracted from the heatmap
![Population number extracted from the heatmap](https://github.com/mungpeter/GAMD_scripts/blob/master/Examples/c_orig_result/extracted_population.png)

Double-mutant thermodynamics cycle
![Double mutant thermodynamics cycle using the population data](https://github.com/mungpeter/GAMD_scripts/blob/master/Examples/c_orig_result/double_mutant_thermo_cycle.png)

##############
- **Running scripts with calculations that are used by the main scripts**
```
gamd_reweight-1d.py   # from original script
gamd_reweight_2d2.py  # modified from gamd_reweight-2d.py, settings are handled by 2_gamd_reweight_run2d.csh
```
- These are the acutal analysis scripts that **2_gamd_reweight_runX.py** scripts call. They are now modified for use with _Python 3.6+_.
- They are modernized from the original analysis script (python2.6-, ca.2014) to enhance the image generation capacity.

```
gamd_reweight-2d.py   # from original script
```
- _gamd_reweight-2d.py_ is the precursor script and is not used by 2_gamd_reweight_run2d.py

###################################################################################
- **Required packages:**
```

anmber        # 16+
cpptraj       # Github version

csh/tcsh      # shell

python        # stable: 3.7.2
  numpy       # stable: 1.16.4
  scipy       # stable: 1.2.1
  pandas      # stable: 0.24.1+
  csv         # stable: 1.0
  matplotlib  # stable: 3.1.0
  tqdm        # stable: 4.32.1
  pathos      # stable: 0.2.3
  argparse    # stable: 1.1 
```
