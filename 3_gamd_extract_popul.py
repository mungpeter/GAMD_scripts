#!/usr/bin/env python3

import sys,os,re
import glob
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathos import multiprocessing
from argparse import ArgumentParser

if int(pd.__version__.split('.')[1]) < 24:
  sys.exit('\n  ERROR: Require Pandas 0.24+ \n')

##########################################################################
#
#  Peter M.U. Ung @ Yale/MSSM
#
#  v.1  19.04.03
#  v.2  19.04.25  do multiple collections with list of X/Y ranges
#  v.3  19.12.09  pandas vectorization in counting
#
#  print out the % population in trajectory in the ranges defined by X-Y axis
#
##  where the data is in 2D-matrix, and bin is defined by the upper end of
##  the range, i.e. range is [0,1,2,3] but data is binned in [0-1, 1-2, 2-3]
#
#####
# example input for -list file: 
  # 0 1 , 3 5.5         # population_1: x-range [0, 1] ; y-range [3, 5.5]
  # 3.5 6.0 ; 3 5.5     # population_2: x-range [3.5,6.0]; y-range [3, 5.5]
  # 3.5 5.0 : 10 15     # population_2: x-range [3.5,5.0]; y-range [10, 15]
#####  
#
#  e.g.> x.py -in     rmsd_file.txt distance_file.txt
#             -col    2
#             -list   multi_ranges.list
#               *or*
#             -xrange 1.0 3.0
#             -yrange 0.0 4.0
#
##########################################################################

def main():

  args = cmdlineparse()     # command line input

  if args.col is None:
    args.col = 2
  if args.skip is None:     # skip X no. of line, if not defined
    args.skip = 0
  if args.pwd is None:
    pwd = os.getcwd()
  else:
    pwd = args.pwd
  print('\033[34m## Working Directory ##\033[0m\n{0}\n'.format(pwd))   # working directory

######################

  # Parse input data then count Frequency of Occurrence in the bins
  Data      = LoadFileToArray(args.input,args.col, skip=int(args.skip),pwd=pwd)
  Frequency = CollectPopulation(data=Data)

  ## Print out results
  if args.xylist:   # manual selection of data
    lpwd = glob.glob(pwd+'/'+args.xylist)[0]
    print("\033[34m## List of data selection ##\033[0m\n{}\n".format(lpwd))

    with open(lpwd, 'r') as fi:
      Ranges = [[x.split() for x in re.split(',|;|:',l.strip())] for l in fi
                              if not re.search('^#', l) ]

    mpi = multiprocessing.Pool(processes=len(Ranges))
    pf  = mpi.map(Frequency, Ranges)
    mpi.close()
    mpi.join()
#    pf = [Frequency(x) for x in Ranges]
  else:
    print('x-range: ', args.x_range)
    print('y-range: ', args.y_range)

    pf = [Frequency( [args.x_range, args.y_range] )]
  
  sum = 0
  for l in pf:
    print(l[0])
    sum += l[1]
  print('# \033[4m% Total Population:\033[0m {:.2f}\n'.format(sum))


##########################################################################
## Count the population found in the specified X/Y ranges to count occurrence
class CollectPopulation(object):
  def __init__( self, data=None ):
    self.data = data

  def __call__( self, ranges ):
    return self.count_in_ranges( ranges )

#########
  def count_in_ranges( self, ranges ):
    x_range, y_range = ranges

    if not x_range or not y_range:
      sys.exit('  ## ERROR: One of the X/Y-ranges is not specified')
    else:
      Xrange = np.array(x_range, dtype=np.float32)
      Yrange = np.array(y_range, dtype=np.float32)

    ## Pandas vectorization, ~100ms for 320K items
    ## create a T/F list of occurrnace in the given X/Y-range as (within)
    ## then create a dataframe with all zeros as long as dataframe
    ## replace the zeros if corresponding position in (within) is True
    df = pd.DataFrame(self.data, columns=['x', 'y'])
    within = ((df.x.to_numpy() >= Xrange[0]) & (df.x.to_numpy() < Xrange[-1]) &
              (df.y.to_numpy() >= Yrange[0]) & (df.y.to_numpy() < Yrange[-1]) )
    occ_df = pd.DataFrame( {'0': [0]*len(df.index)} )
    occ_df[within == True] = 1
    count  = occ_df['0'].sum()
    
    ## Old method: typical loop of loop: ~5x slower than vectorization
#    count = 0
#    for pt in self.data:
#      if pt[0] >= Xrange[0] and pt[0] < Xrange[-1]:
#        if pt[1] >= Yrange[0] and pt[1] < Yrange[-1]:
#          count = count + 1

    percent = count*100/len(df.index)
    pf = '#\033[34m X-range\033[0m [ {0:.1f} , {1:.1f} ] | \033[34mY-range\033[0m [ {2:.1f} , {3:.1f} ]\n \033[4m% Population:\033[0m {4:.2f}\n'.format(
      Xrange[0], Xrange[1], Yrange[0], Yrange[1], percent )

    return pf, percent


##########################################################################
# Load data files with numpy, automatically skip comment lines with '#', 
# convert input to array of columns; use specified column as input 
#skip = 1    # in case too much data, slice the data array by this number
def LoadFileToArray( files, col, pwd='', skip=0 ):

    obj = LoadData(col=col, skip=skip, pwd=pwd)
    mpi = multiprocessing.Pool(processes = 2)
    Tmp = mpi.map(obj, files)
    mpi.close()
    mpi.join()
#    Tmp = [obj(x) for x in files]   

    if len(Tmp[0]) != len(Tmp[1]):
        sys.exit('## ERROR: files do not match in length:\n{0}\t{1}\n{2}\t{3}'.format(
                  files[0], len(Tmp[0]), files[1], len(Tmp[1])))

    loaded = np.array(list(zip(Tmp[0], Tmp[1])))
    return loaded


##########################################################################
class LoadData(object):
  def __init__( self, col='', pwd='', skip=0 ):
    self.col  = col
    self.pwd  = pwd
    self.skip = skip+1
  
  def __call__( self, infile ):
    return self.load_txt( infile )
  
##########
  # default 'comment' of np.loadtxt = '#'
  def load_txt( self, infile ):
    print(self.pwd+'/'+infile)
    fpwd = glob.glob(self.pwd+'/'+infile)[0]
    load = np.loadtxt(fpwd, comments='#', usecols=[int(self.col)-1] )[::self.skip]
    print("\033[34m## DATA LOADED:\033[0m {0} lines\n> {1}\n".format(len(load), fpwd))
    return load


##########################################################################
def cmdlineparse():
  p = ArgumentParser(description="command line arguments")
  p.add_argument("-in", dest="input", required=True, nargs="+",
                 help="Input data files", metavar="<X-input Y-input>")

  p.add_argument('-col', dest='col', required=False,
                 help='Column in file to read (def: 2)', metavar='<column>')

  p.add_argument("-pwd", dest='pwd', required=False,
                 help="Path of working directory", metavar='<work directory>')
  p.add_argument("-skip", dest='skip', required=False,
                 help="Skip every other X input lines", metavar='<skip line>')

  p.add_argument("-xrange", dest="x_range", required=False, nargs='+',
                 help="Range of X-dimension to collect", metavar="<xmin xmax>")
  p.add_argument("-yrange", dest="y_range", required=False, nargs='+',
                 help="Range of Y-dimension to collect", metavar="<ymin ymax>")
  p.add_argument("-list", dest="xylist", required=False,
                 help="List of X- and Y-ranges to collect, x-/y-ranges separated by delimiters ',|;|:'", metavar="<range list>")

  #####
  # example input for -list file: 
  # 0 1 , 3 5.5         # population_1: x-range [0, 1] ; y-range [3, 5.5]
  # 3.5 6.0 ; 3 5.5     # population_2: x-range [3.5,6.0]; y-range [3, 5.5]
  # 3.5 5.0 : 10 15     # population_2: x-range [3.5,5.0]; y-range [10, 15]
  #####  

  args = p.parse_args()
  return args

##########################################################################
if __name__ == '__main__':
  main()
