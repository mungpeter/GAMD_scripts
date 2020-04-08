#!/bin/csh -f

if ($#argv != 4) then
  echo
  echo "   ## Usage: x.csh [LSF template] [.rst] [Repeat No.] [Curr No.] ##"
  echo
  echo "        e.g: x.csh hpc.lsf 3E_equilibrate.1.rst 10 6"
  exit
endif

set lsf    = $argv[1]	# Template AMBER simulation script
set prmcrd = $argv[2]	# MD restart file
set repeat = $argv[3]	# Total number of repeat
set curRun = $argv[4]	# Current iteration

if ($curRun <= $repeat) then

  sed "s/REFFILE/$lsf/" $lsf | \
  sed "s/PRMCRD/$prmcrd/"    | \
  sed "s/REPEAT/$repeat/"    | \
  sed "s/CURRUN/$curRun/"      \
       > $lsf.$curRun
  bsub < $lsf.$curRun

endif
exit
