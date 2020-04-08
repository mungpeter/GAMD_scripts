#!/usr/bin/env python3

############################################################################

# Peter M.U. Ung @ Yale

# v1.
# v2.   19.10.23 - cleanup

#  main script to compile data generated from GAMD simulations. 
#  take in the GAMD-generated energy file and proteinstructural metrics 
#  as input to generate 2D heatmap of potential-of-mean-force free energy
#  landscape and then 2D conformational population landscape images.
#
#  require at least 30,000 data points from 300ns simulations to have enough
#  sampling to generate good results
#
#  image parameters can be fine-tuned

############################################################################

import sys,os
import re,glob
import subprocess
import numpy as np
from tqdm import tqdm
from pathos import multiprocessing
from argparse import ArgumentParser
from gamd_reweight_2d2 import gamd_reweight_2d

def main():

  args = cmdparse()

  if args.pwd is None:
    pwd = os.getcwd()
  else:
    pwd = args.pwd
  print('\033[34m## Working Directory ##\033[0m\n{0}\n'.format(pwd))   # working directory
  if args.dir:
    pwd = pwd+'/'+args.dir

  if args.job is None:
    job = 'amdweight_CE'
  else:
    job = args.job

  if args.col is None:
    args.col = 2

  # X-/Y-axes bin sizes are generally controlled by args.bin, but if specific
  # flag is used, each axis' bin size is changed accordingly
  if not args.bin:
    xbin, ybin = '0.5', '0.5'
  else:
    xbin, ybin = args.bin, args.bin
  if args.Xbin:
    xbin = args.Xbin
  if args.Ybin:
    ybin = args.Ybin
  binz = [xbin, ybin]

##  SET TEMPERATURE
  if args.T:
    T = float(args.T)
  else:
    T = 310   # simulation temperature

##  SET MAX ENERGY FOR ALL INFINITY VALUES
  if args.Emax:
    e_max = float(str(args.Emax))
  else:
    e_max = 4.

##  SET MAX POPULATION
  if args.Pmax:
    p_max = float(str(args.Pmax))
  else:
    p_max = 5.

##  SET HISTOGRAM CUTOFF
  if args.cut:
    cut = int(args.cut)
  else:
    cut = 10  # minimum number of configurations in one bin

##  SET ORDER of McLaurin series expansion
  if args.order:
    order = int(args.order)
  else :
    order = 10	# default

##  SET flag for Gaussian fitting of deltaV
  if args.fit:
    fit = args.fit
  else :
    fit = False	# simulation temperature

## SET FIG CONTOUR SMOOTHENING LEVEL
  if args.smooth:
    smo = float(args.smooth)
  else:
    smo = None

## SET FIGURE RESOLUTION DPI
  if args.dpi:
    dpi = int(args.dpi)
  else:
    dpi = 200

## SET FIGURE FORMAT
  if args.img:
    img = args.img
  else:
    img = 'png'

## SET FIG CONTOUR NUMBER per Emax integer
## default is 4 when args.bin is 0.5
  if args.c_step:
    c_step = float(args.c_step)
  else:
    c_step = 4.


############################
  ## read in input files and x/y axis names
  with open(args.input, 'r') as fi:
    inplist = [ l.rstrip().split() for l in fi 
                  if l.rstrip() and not re.search('^#',l) ]

  ## read in the each input file once, put them into dictionary
  inpfiles = set([item for sublist in list(zip(*inplist))[:2] for item in sublist ])

  Rf = LoadData(pwd=pwd, col=args.col, skip=1)
  mpi = multiprocessing.Pool()
#  InpFiles = [x for x in tqdm(mpi.imap_unordered(
#                         Rf, inpfiles), total=len(inpfiles))]
  InpFiles = [Rf(x) for x in inpfiles]
  Dict = {}
  for f in InpFiles:
    Dict[list(f.keys())[0]] = f[list(f.keys())[0]]
  
  weights, dV = WeightParse( args.gamd, job, pwd=pwd, skip=1)

  gmd = GAMD( gamd=weights, dV=dV, data=Dict, 
              job=args.job, e_max=e_max, p_max=p_max, dimn='2D',
              binz=binz, T=T, xdim=args.Xdim, ydim=args.Ydim, 
              cut=cut, order=order, fit=fit,
              c_step=c_step, smooth=smo, dpi=dpi, img=img )

#  # cumulant expansion #
#  m = [gmd.cumulant_expansion(x) for x in tqdm(inplist)]
  m = [x for x in tqdm( mpi.imap_unordered(gmd.cumulant_expansion, inplist),
                                  total=len(inplist) )]
#  # maclaurin series #
#  n = [x for x in tqdm( mpi.imap(gmd.maclaurin_series, inplist), 
#                                 total=len(inplist) )]
#  # reweighting #
#  o = [x for x in tqdm( mpi.imap(gmd.reweighting, inplist),
#                                 total=len(inplist) )]

#  #no weighting #
#  p = [x for x in tqdm( mpi.imap(gmd.no_weighting, inplist),
#                                 total=len(inplist) )]

#  gmd.amd_time(inplist[0])


##########################################################################
class GAMD(object):
  def __init__( self, gamd=[], dV=[], data={}, job='', e_max='', p_max='',
                      binz=[], T='', dimn='', xdim=[], ydim=[], 
                      cut='', order='', fit='', 
                      c_step='', dpi='', img='', smooth='' ):
    self.gamd = gamd        # GAMD energy result
    self.dV   = dV          # GAMD dV from read in
    self.data = data        # dictionary of input data

    self.binz = binz        # general bin size to cluster
    self.T    = T           # temperature

    self.p_max= p_max       # max population percentage cutoff
    self.dimn = dimn        # dimension of input data (1D or 2D)
    self.xdim = xdim        # user-specified x-axis scale
    self.ydim = ydim        # user-specified y-axis scale
    self.job  = job         # GAMD analysis type
    self.e_max= e_max       # Max Energy cutoff in graph
    self.cut  = cut         # cutoff for histogram
    self.order= order       # McLaurin series expansion order
    self.fit  = fit         # flag for Gaussian fitting of deltaV
    self.smo  = smooth      # smoothen contour data with spline
    self.dpi  = dpi         # figure DPI
    self.img  = img         # image format
    self.c_step = c_step    # number of contour level between integer number
    
################
  def _check_name(self, rawfile):
    tmp = rawfile   # filename1, filename2, x-label, y-label

    if len(tmp) == 4: # case: 2 file input, x-lab, y-lab
      name = tmp[0].split('.txt')
      part = tmp[1].split('.txt')[0].split('.')[-1]
      infi = name[0]+'.'+part  # unique name by combining 2 different files
      return [tmp[0]+' '+tmp[1], infi, '.txt'+name[1], tmp[-2], tmp[-1]]

    if len(tmp) == 3: # case: 1 file input, x-lab, y-lab
      name =tmp[0].split('.txt')
      return [tmp[0], name, '.txt'+name[1], tmp[-2], tmp[-1]]

    if len(tmp) == 2: # case: 1 file input, x-lab
      name =tmp[0].split('.txt')
      return [tmp[0], name[0], '.txt'+name[1], tmp[-1], '' ]

    else:   # case: 1 file input
      name = tmp.split('.txt')
      return [tmp, name[0], '.txt'+name[1], '', '']

################
  def cumulant_expansion( self, rawfile ):
    in_file, combin, subfx, xlab, ylab = self._check_name(rawfile)
    inp = in_file.split()
    job = 'amdweight_CE'

    if len(self.data[inp[0]]) != len(self.data[inp[1]]):
      sys.exit('\033[31m## ERROR: files do not match in length:\033[0m\n{0}\t{1}\n{2}\t{3}'.format(
                inp[0], len(self.data[inp[0]]), inp[1], len(self.data[inp[1]])))
    elif len(self.data[inp[0]]) != len(self.gamd):
      sys.exit('\033[31m## ERROR: Input and GAMD Weight files do not match in length:\033[0m\n{0}\t{1}\n{2}\t{3}'.format(
                inp[0], len(self.data[inp[0]]), 'GAMD weights', len(self.gamd)) )
    else:
      data_2d = np.array(list(zip( self.data[inp[0]], self.data[inp[1]] )))


    gamd_reweight_2d( combiname=combin, data=data_2d, weights=self.gamd, dV=self.dV,
                      job=job, T=self.T, Xdim=self.xdim, Ydim=self.ydim, 
                      binz=self.binz, e_max=self.e_max, p_max=self.p_max,
                      cutoff=self.cut, 
                      order=self.order, fit=self.fit, Xlab=xlab, Ylab=ylab, 
                      dpi=self.dpi, img=self.img, c_step=self.c_step, smooth=self.smo )

#    os.system('python {} -input {} -col {} -weight {} -job {} {} {} {} {} {} {} {} {} {} {} {} {} | tee -a {}.{}-reweight_CE.log'.format( 
#              self.pyrw, in_file, self.col, self.gamd, job, self.binz, self.T,
#              self.xdim, self.ydim, self.cut, self.e_max, self.p_max, 
#              xlab, ylab, self.dpi, self.img,
#              self.c_step, self.smo, self.pwd,    name, self.dimn))  

    os.system('mv pmf-c1-{0}.xvg pmf-{0}.{1}-CE1.xvg'.format(
              combin, self.dimn))
    os.system('mv pmf-c2-{0}.xvg pmf-{0}.{1}-CE2.xvg'.format(
              combin, self.dimn))
    os.system('mv pmf-c3-{0}.xvg pmf-{0}.{1}-CE3.xvg'.format(
              combin, self.dimn))


################
  def maclaurin_series(self, rawfile):
    in_file, combin, subfx, xlab, ylab = self._check_name(rawfile)
    inp = in_file.split()
    job='amdweight_MC'

    if len(self.data[inp[0]]) != len(self.data[inp[1]]):
      sys.exit('## ERROR: files do not match in length:\n{0}\t{1}\n{2}\t{3}'.format(
                inp[0], len(self.data[inp[0]]), inp[1], len(self.data[inp[1]])))
    elif len(self.data[inp[0]]) != len(self.gamd):
      sys.exit('## ERROR: Input and GAMD Weight files do not match in length:\n{0}\t{1}\n{2}\t{3}'.format(
                inp[0], len(self.data[inp[0]]), 'GAMD weights', len(self.gamd)) )
    else:
      data_2d = np.array(list(zip( self.data[inp[0]], self.data[inp[1]] )))

    gamd_reweight_2d( combiname=combin, data=data_2d, gamd=self.gamd, dV=self.dV,
                      job=job, T=self.T, Xdim=self.xdim, Ydim=self.ydim, 
                      binz=self.binz, e_max=self.e_max, p_max=self.p_max,
                      cutoff=self.cut, 
                      order=self.order, fit=self.fit, Xlab=xlab, Ylab=ylab, 
                      dpi=self.dpi, img=self.img, c_step=self.c_step, smooth=self.smo )

#    os.system('python {} -input {} -col {} -weight {} -job {} {} {} -order 10 {} {} {} {} {} {} {} {} {} {} | tee -a {}.{}-reweight_MC.log'.format( 
#              self.pyrw, in_file, self.col, self.gamd, job, self.binz, self.T,
#              self.xdim, self.ydim, self.cut, self.e_max, self.p_max, 
#              xlab, ylab, self.dpi, self.img,
#              self.c_step, self.smo, self.pwd,   name, self.dimn))

    os.system('mv pmf-{0}.xvg pmf-{0}.{1}-MCorder10.xvg'.format(
              combin, self.dimn))


################
  def reweighting(self, rawfile):
    in_file, combin, subfx, xlab, ylab = self._check_name(rawfile)
    inp = in_file.split()
    job='amdweight'

    if len(self.data[inp[0]]) != len(self.data[inp[1]]):
      sys.exit('## ERROR: files do not match in length:\n{0}\t{1}\n{2}\t{3}'.format(
                inp[0], len(self.data[inp[0]]), inp[1], len(self.data[inp[1]])))
    elif len(self.data[inp[0]]) != len(self.gamd):
      sys.exit('## ERROR: Input and GAMD Weight files do not match in length:\n{0}\t{1}\n{2}\t{3}'.format(
                inp[0], len(self.data[inp[0]]), 'GAMD weights', len(self.gamd)) )
    else:
      data_2d = np.array(list(zip( self.data[inp[0]], self.data[inp[1]] )))

    gamd_reweight_2d( combiname=combin, data=data_2d, gamd=self.gamd, dV=self.dV,
                      job=job, T=self.T, Xdim=self.xdim, Ydim=self.ydim, 
                      binz=self.binz, e_max=self.e_max, p_max=self.p_max,
                      cutoff=self.cut, 
                      order=self.order, fit=self.fit, Xlab=xlab, Ylab=ylab, 
                      dpi=self.dpi, img=self.img, c_step=self.c_step, smooth=self.smo )

#    os.system('python {} -input {} -col {} -weight {} -job {} {} {} {} {} {} {} {} {} {} {} {} {} | tee -a {}.{}-reweight.log'.format( 
#              self.pyrw, in_file, self.col, self.gamd, job, self.binz, self.T,
#              self.xdim, self.ydim, self.cut, self.e_max, self.p_max,
#              xlab, ylab, self.dpi, self.img,
#              self.c_step, self.smo, self.pwd,    name, self.dimn))

    os.system('mv pmf-{0}.xvg pmf-{0}.{1}-reweight.xvg'.format(
              combin, self.dimn))


################
  def no_weighting(self, rawfile):
    in_file, combin, subfx, xlab, ylab = self._check_name(rawfile)
    inp = in_file.split()
    job='noweight'

    if len(self.data[inp[0]]) != len(self.data[inp[1]]):
      sys.exit('## ERROR: files do not match in length:\n{0}\t{1}\n{2}\t{3}'.format(
                inp[0], len(self.data[inp[0]]), inp[1], len(self.data[inp[1]])))
    elif len(self.data[inp[0]]) != len(self.gamd):
      sys.exit('## ERROR: Input and GAMD Weight files do not match in length:\n{0}\t{1}\n{2}\t{3}'.format(
                inp[0], len(self.data[inp[0]]), 'GAMD weights', len(self.gamd)) )
    else:
      data_2d = np.array(list(zip( self.data[inp[0]], self.data[inp[1]] )))

    gamd_reweight_2d( combiname=combin, data=data_2d, gamd=self.gamd, dV=self.dV,
                      job=job, T=self.T, Xdim=self.xdim, Ydim=self.ydim, 
                      binz=self.binz, e_max=self.e_max, p_max=self.p_max,
                      cutoff=self.cut, 
                      order=self.order, fit=self.fit, Xlab=xlab, Ylab=ylab, 
                      dpi=self.dpi, img=self.img, c_step=self.c_step, smooth=self.smo )

#    os.system('python {} -input {} -col {} -weight {} -job {} {} {} {} {} {} {} {} {} {} {} {} {} | tee -a {}.{}-no_weight.log'.format( 
#              self.pyrw, in_file, self.col, self.gamd, job, self.binz, self.T,
#              self.xdim, self.ydim, self.cut, self.e_max, self.p_max,
#              xlab, ylab, self.dpi, self.img,
#              self.c_step, self.smo, self.pwd,   name, self.dimn))

    os.system('mv pmf-{0}.xvg pmf-{0}.{1}-noweight.xvg'.format(
              combin, self.dimn))


################
  def amd_time(self, rawfile):
    in_file, combin, subfx, xlab, ylab = self._check_name(rawfile)
    inp = in_file.split()
    job='amd_time'

    if len(self.data[inp[0]]) != len(self.data[inp[1]]):
      sys.exit('## ERROR: files do not match in length:\n{0}\t{1}\n{2}\t{3}'.format(
                inp[0], len(self.data[inp[0]]), inp[1], len(self.data[inp[1]])))
    elif len(self.data[inp[0]]) != len(self.gamd):
      sys.exit('## ERROR: Input and GAMD Weight files do not match in length:\n{0}\t{1}\n{2}\t{3}'.format(
                inp[0], len(self.data[inp[0]]), 'GAMD weights', len(self.gamd)) )
    else:
      data_2d = np.array(list(zip( self.data[inp[0]], self.data[inp[1]] )))

    gamd_reweight_2d( combiname=combin, data=data_2d, gamd=self.gamd, dV=self.dV,
                      job=job, T=self.T, Xdim=self.xdim, Ydim=self.ydim, 
                      binz=self.binz, e_max=self.e_max, p_max=self.p_max,
                      cutoff=self.cut, 
                      order=self.order, fit=self.fit, Xlab=xlab, Ylab=ylab, 
                      dpi=self.dpi, img=self.img, c_step=self.c_step, smooth=self.smo )

#    os.system('python {} -input {} -col {} -weight {} -job {} {} {} | tee -a {}.{}-amd_time.log'.format( 
#              self.pyrw, in_file, self.gamd, job, self.binz, self.T, self.pwd,
#              name, self.dimn))


##########################################################################
# Load data files with numpy, automatically skip comment lines with '#', 
# convert input to array of columns; use specified column as input 
#skip = 1    # in case too much data, slice the data array by this number

class LoadData(object):
  def __init__( self, pwd='', col='', skip='' ):
    self.pwd  = pwd
    self.col  = col
    self.skip = skip

  def __call__( self, infile ):
    return self.load_txt( infile )

  # default 'comment' of np.loadtxt = '#'
  def load_txt( self, infile ):
    print('\033[34m## Looking into directory for data: \033[0m')
    print(glob.glob(self.pwd))
    print(infile)
    try:
      print(glob.glob(self.pwd+'/'+infile))
    except IndexError:
      sys.exit('\033[31m  ERROR: Either filename is incorrect, or path has special characters\033[30m')  

    fpwd = glob.glob(self.pwd+'/'+infile)[0]
    load = np.loadtxt(fpwd, comments='#',
                      usecols=[int(self.col)-1])[::self.skip]
    print("\033[34mDATA LOADED:\033[0m {0} lines\n".format(len(load)))
    return { infile: load }

##########################################################################
def WeightParse( gamd_file, job, pwd='', skip=1 ):
  print('\033[34m## Looking into directory for GAMD Energy file: \033[0m')
  print('\033[34m   Using Reweighting method: \033[0m'+job)
  if job == "weighthist":
    try:
      gpwd = glob.glob(pwd+'/'+gamd_file)[0]
    except IndexError:
      print(pwd+'/'+gamd_file)
      sys.exit('\033[31m  ERROR: Either filename is incorrect, or path has special characters\033[30m')
    print(gpwd)
    data = np.loadtxt(gpwd, comments='#')[::skip]
    print("\033[34mDATA LOADED:\033[0m {0} lines\n".format(len(data)))
    weights = data[:,0]
    dV   = np.zeros(len(weights))
  elif job == "amd_time" or job == "amd_dV" or job == "amdweight" or job == "amdweight_MC" or job == "amdweight_CE" :
    try:
      gpwd = glob.glob(pwd+'/'+gamd_file)[0]
    except IndexError:
      print(pwd+'/'+gamd_file)
      sys.exit('\033[31m  ERROR: Either filename is incorrect, or path has special characters\033[30m')
    print(gpwd)
    data = np.loadtxt(gpwd, comments='#')[::skip]
    print("\033[34mDATA LOADED:\033[0m {0} lines\n".format(len(data)))
    weights = np.exp(data[:,0])
    dV = data[:,2]
  elif job == "noweight":
    data = np.loadtxt(gamd_file, comments='#')[::skip]
    weights = np.zeros(len(data))
    weights = weights + 1
    dV      = np.zeros(len(data))
  else:
    print( "ERROR JOBTYPE"+ job + " NOT RECOGNIZED")
    del data
    del weights
    del dV
  return weights, dV


###########################################################################

def cmdparse():
  parser = ArgumentParser(description="command line arguments")

  parser.add_argument("-inp", dest="input", required=True, 
                      help="list of input files (format: x-axis-file y-axis-file x-axis-name y-axis-name)", metavar="<input list>")
  parser.add_argument("-gamd", dest="gamd", required=True, 
                      help="GAMD weight file (no path needed)", metavar="<GAMD weight file>")

  parser.add_argument('-col', dest='col', required=False,
                      help='column to be read (def: 2)', metavar='<column number>')

  parser.add_argument('-pwd', dest='pwd', required=False,
                      help="Path of working directory, beware of characters' \()' (def: ./)", metavar='<work directory>')
  parser.add_argument('-dir', dest='dir', required=False,
                      help='path to input files only, on top of -pwd (def: ./)', metavar='<path>')

  parser.add_argument('-job', dest='job', required=False,
                      help='Reweighting method: <noweight>, <weighthist>, <amd_time>, <amd_dV>, <amdweight>, <amdweight_MC>, <amdweight_CE>: (Default: amdweight_CE)', metavar='<Job type reweighting method>')
  parser.add_argument("-bin", dest="bin", required=False,
                      help="General bin size (def: 0.5)", metavar="<bin size>")
  parser.add_argument("-Xbin", dest="Xbin", required=False,
                      help="Bin size on X-axis, supersede -bin (def: auto)", metavar="<Xbin size>")
  parser.add_argument("-Ybin", dest="Ybin", required=False,
                      help="Bin size on Y-axis, supersede -bin (def: auto)", metavar="<bin size>")
  parser.add_argument("-Emax", dest="Emax", required=False,
                      help="Max energy above 0 kcal/mol (def: 4)", metavar="<Emax>")

  parser.add_argument("-Pmax", dest="Pmax", required=False,
                      help="Max Population percentage (def: 5)", metavar="<Pmax>")

  parser.add_argument("-temp", dest="T", required=False, 
                      help="Temperature K (def: 310)",metavar="<temperature>")
  parser.add_argument("-Xdim", dest="Xdim", required=False, nargs="+", 
                      help="Data range of X-dimension (def: auto)", metavar="<Xmin Xmax>")
  parser.add_argument("-Ydim", dest="Ydim", required=False, nargs="+", 
                      help="Data range of Y-dimension (def: auto)", metavar="<Ymin Ymax>")
  parser.add_argument("-histcut", dest="cut", required=False,
                      help="Data cutoff (def: 25)", metavar="<histo cutoff>")
  parser.add_argument("-order", dest="order", required=False,
                      help="Order of MC formula (def: 10)", metavar="<Mc order>")
  parser.add_argument("-fit", dest="fit", required=False,
                      help="Fitting data (def: False)", metavar="<Fit switch>")


  parser.add_argument('-contour', dest='c_step', required=False,
                      help='Figure Contour per integer (def: 4)', metavar='<contour>')
  parser.add_argument('-smooth', dest='smooth', required=False,
                      help='Contour smoothening (def: None | Sug: 1.15)', metavar='<smooth>')
  parser.add_argument('-dpi', dest='dpi', required=False,
                      help='Figure Resolution (def: 200)', metavar='<dpi>')
  parser.add_argument('-img', dest='img', required=False,
                      help='Figure image format: png|svg|eps|pdf (def: png)', metavar='<img>')

  args = parser.parse_args()
  return(args)


##########################################################################

if __name__ == '__main__':
  main()

