#!/usr/bin/env python3

###########################################################################
#
#  Peter M.U. Ung @ MSSM/Yale
#
#  v1.0  19.05.??
#  v2.0  19.06.21  convert to py3 syntax, encode/decode byte strings
#
#  function: prepare GAMD.log files from multiple sources
#
#            This script will take in all specified gamd.log.bz2 files in
#            main-folders/sub-folders and comply a summarized gamd.log.bz2
#            output, while removing the commented lines and the first line
#            of trajectory with 0000000.
#            Opted to use manual IO instead of numpy getfromtxt since it is
#            much faster to do it without numpy.
#
# e.g.>  x.py \
#            -mf 1_v16 2_v16a 3_v16r 4_v16ar \
#            -sf 1_run 2_run                 \
#            -log gamd.1.log.bz2 gamd.2.log.bz2 gamd.3.log.bz2 \
#            -o fgf21_variants_gamd          \
#            -t 310 -s 1
#
# ** directory structure
# ----1_v16____1_run_____gamd.1.log.bz2
#   |    |         |_____gamd.2.log.bz2
#   |    |         |_____gamd.3.log.bz2
#   |    |
#   |    |______2_run____gamd.1.log.bz2
#   |               |____ ...
#   |
#   |-2_v16a___1_run_____ ...
#   |    |_____2_run_____ ...
#   |
#   |-3_v16r___ ...
#   |-4_v16ar__ ...
#
###########################################################################
import sys,re
import gzip,bz2
import numpy as np
from argparse import ArgumentParser


def main():

  Engs = []
  args = cmdparse()

  if not args.step:
    step = 1
  else:
    step = int(args.step)

  if not args.T:
    T = 310
  else:
    T = float(args.T)

####################################
  num = 0
  for folder in args.folder:
    for sub in args.sub:
      for log in args.log:
        Engs = Engs + process_gamd('{}/{}/{}'.format(folder,sub,log),step,T)
        num = num+1

######################################
## Output, always gives .bz2 format result
## py3.5+ requires encoding strings to byte for bz2 output
  num = 0
  outfile = args.out+'.dat.bz2'

  with bz2.BZ2File(outfile, 'wb', compresslevel=9) as fo:
    for l in Engs:
      # column 2 (timestamp) is not important
      fo.write('{0:.8f}\t{1}\t{2:.8f}\n'.format(l[0], l[1], l[2]).encode())
      num = num+1
    print('\033[34mData point read: \033[0m'+str(num))
  
  
##########################################################################
## py3.5+ requires decoding the byte to string
def process_gamd(infile, step, T):
  
  data = []
  with file_handle(infile) as fi:
    print('> '+infile)

    for m in fi:
      if re.search(rb'#', m):
        continue
      l = m.decode('utf-8').split()
      if int(float(l[2])) == 0 and int(float(l[3])) == 0:
        print('\033[34m  > remove 0.00000\033[0m')
      else:
        boost_en = float(l[7])+float(l[6])
        reweight = boost_en/(0.0019872036*T)
        data.append([reweight, int(l[1]), boost_en])

  Part = data[::step]
  del(data)
  return(Part)

##########################################################################
## numpy readin method but twice as slow as normal way - won't be used
def process_gamd_np(infile, step, T):

  data = []
  for l in np.loadtxt(infile):
    if int(l[2]) == 0 and int(l[3]) == 0:
      print(infile)
    else:
      boost_en = l[7]+l[6]
      reweight = boost_en/(0.0019872036*T)
      data.append([reweight, int(l[1]), boost_en])

  Part = data[::step]
  del(data)
  return(Part)

##########################################################################
def file_handle(file_name):
  if re.search(r'.gz$', file_name):
    handle = gzip.open(file_name, 'r')
  elif re.search(r'.bz2$', file_name):
    handle = bz2.BZ2File(file_name, 'r')
  else:
    handle = (file_name)
  return handle


##########################################################################
def cmdparse():
  p = ArgumentParser(description='## prepare GAMD.log files from multiple sources\n')
  
  p.add_argument('-mf', dest='folder', required=True, nargs='+',
                  help='main folders (in order)', metavar='<fold1 fold2 ...>')
  p.add_argument('-sf', dest='sub', required=True, nargs='+',
                  help='sub folders (in order)', metavar='<sub1 sub2 ...>')
  p.add_argument('-log', dest='log', required=True, nargs='+',
                  help='gamd logs (in order)', metavar='<log1 log2 ...>')
  p.add_argument('-o', dest='out', required=True,
                  help='output file prefix', metavar='<filename>')

  p.add_argument('-t', dest='T', required=False,
                  help='temperature K (def: 310)', metavar='<Temperature>')
  p.add_argument('-s', dest='step', required=False,
                  help='Step interval (def: 0)', metavar='<stepsize>')

  args = p.parse_args()
  return(args)


##########################################################################

if __name__ == '__main__':
  main()
