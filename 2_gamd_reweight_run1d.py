#!/usr/bin/env python3

import sys,os
import re,glob
import subprocess
from tqdm import tqdm
from pathos import multiprocessing
from argparse import ArgumentParser

gamd1d = '/home/pmung/Dropbox/9_scripts/2_MD/1_amber/4_gamd/gamd_reweight-1d.py'
#gamd1d = './gamd_reweight-1d.py'
#gamd2d = '/home/pmung/Dropbox/9_scripts/2_MD/1_amber/4_gamd/gamd_reweight-2d.py'
#gamd2d = './gamd_reweight-2d.py'

def main():

  dimn = '1D'
  program = gamd1d

############
  args = cmdparse()

  if not args.pwd:
    pwd = '-pwd .'
  else:
    pwd = '-pwd '+args.pwd
  if args.dir:
    pwd = pwd+'/'+args.dir

  if args.col is None:
    args.col = '2'

  # X-/Y-axes bin sizes are generally controlled by args.bin, but if specific
  # flag is used, each axis' bin size is changed accordingly
  if not args.bin:
    binz = '-bin 0.5'
  else:
    binz = '-bin '+args.bin

  if args.T:
    T = '-T '+args.T
  else:
    T = '-T 310'

  if args.Emax:
    emax = '-Emax '+args.Emax
  else:
    emax = '-Emax 4.'

  if args.cut:
    cut = '-cutoff '+args.cut
  else:
    cut = ''

  if args.smooth:
    smo = '-smooth '+args.smooth
  else:
    smo = ''

  if args.dpi:
    dpi = '-dpi '+args.dpi
  else:
    dpi = '-dpi 200'

  if args.Xdim:
    tmp  = args.Xdim
    xdim = '-Xdim {0} {1}'.format(tmp[0], tmp[1])
  else:
    xdim = ''

  if args.c_step:
    c_step = '-contour '+args.c_step
  else:
    c_step = '-contour 4'

############################
  mpi = multiprocessing.Pool()
  gmd = GAMD( program=program, gamd=args.gamd, binz=binz, col=args.col,
              pwd=pwd, T=T, dimn=dimn, xdim=xdim, emax=emax, 
              c_step=c_step, smooth=smo, dpi=dpi, cut=cut )

#  m = [gmd.cumulant_expansion(x) for x in open(args.input, 'r')]
  with open(args.input, 'r') as fi:
    inplist = [ l.rstrip() for l in fi 
                  if l.rstrip() and not re.search('^#', l)]
  print(inplist)

  print('\033[34m# cumulant expansion #\033[0m')
#  m = [gmd.cumulant_expansion(x) for x in tqdm(inplist)]
  m = [x for x in tqdm( mpi.imap_unordered(gmd.cumulant_expansion, inplist),
                                 total=len(inplist))]
  print('\033[34m# maclaurin series #\033[0m')
#  n = [x for x in tqdm( mpi.imap(gmd.maclaurin_series, inplist), 
#                                 total=len(inplist))]
  print('\033[34m# reweighting #\033[0m')
#  o = [x for x in tqdm( mpi.imap(gmd.reweighting, inplist),
#                                 total=len(inplist))]

  print('\033[34m# no weighting #\033[0m')
#  p = [x for x in tqdm( mpi.imap(gmd.no_weighting, inplist),
#                                 total=len(inplist))]

#  gmd.amd_time(inplist[0])


##########################################################################
class GAMD(object):
  def __init__(self, program='', gamd='', col='', binz='',
                     pwd='', T='', dimn='', xdim='', emax='', 
                     cut='', c_step='', dpi='', smooth='' ):
    self.gamd = gamd        # GAMD energy result
    self.col  = col         # column of data in input file
    self.binz = binz        # general bin size to cluster
    self.pwd  = pwd         # directory to input data
    self.T    = T           # temperature
    self.dimn = dimn        # dimension of input data (1D or 2D)
    self.xdim = xdim        # user-specified x-axis scale
    self.emax = emax        # Max Energy cutoff in graph
    self.cut  = cut         # cutoff for histogram
    self.smo  = smooth      # smoothen contour data with spline
    self.dpi  = dpi         # figure DPI
    self.pyrw = program     # which program to run, 1D or 2D
    self.c_step = c_step    # number of contour level between integer number
    
################
  def _check_name(self, rawfile):
    inp = rawfile.rstrip() 
    if re.search(r' ', inp):
      tmp = inp.split()	# filename, x-label, y-label

      if len(tmp) == 2: # case: 1 file input, x-lab
        name =tmp[0].split('.txt')
        return [tmp[0], name[0], '.txt'+name[1], '-Xlab '+tmp[-1] ]
      else:
        sys.exit('\033[31m    ERROR: Input list has Maximum 2items per line: [file] [x-label]\033[0m')
    else:   # case: 1 file input
      name = inp.split('.txt')
      return [inp, name[0], '.txt'+name[1], '']

################
  def cumulant_expansion(self, rawfile):
    in_file, name, subfx, xlab = self._check_name(rawfile)
    job = 'amdweight_CE'

#    print('xxxxxxxxxxx', job, self.gamd)
    os.system('python {} -input {} -col {} -gamd {} -job {} {} {} {} {} {} {} {} {} {} {} | tee -a {}.{}-reweight_CE.log'.format( 
              self.pyrw, in_file, self.col, self.gamd, job, self.binz, self.T,
              self.xdim, self.cut, self.emax, xlab, self.dpi,
              self.c_step, self.smo, self.pwd,    name, self.dimn))  

    os.system('mv -v pmf-c1-{0}.xvg pmf-{1}.{2}-CE1.xvg'.format(
              name+subfx, name, self.dimn))
    os.system('mv -v pmf-c2-{0}.xvg pmf-{1}.{2}-CE2.xvg'.format(
              name+subfx, name, self.dimn))
    os.system('mv -v pmf-c3-{0}.xvg pmf-{1}.{2}-CE3.xvg'.format(
              name+subfx, name, self.dimn))
#    os.system('xmgrace -hdevice SVG -hardcopy -nxy pmf-{0}.{1}-CE1.xvg -printfile pmf-{0}.{1}-CE1.svg'.format(name, self.dimn))
    os.system('xmgrace -hdevice SVG -hardcopy -nxy pmf-{0}.{1}-CE2.xvg -printfile pmf-{0}.{1}-CE2.svg'.format(name, self.dimn))
#    os.system('xmgrace -hdevice SVG -hardcopy -nxy pmf-{0}.{1}-CE2.xvg -printfile pmf-{0}.{1}-CE3.svg'.format(name, self.dimn))


################
  def maclaurin_series(self, rawfile):
    in_file, name, subfx, xlab = self._check_name(rawfile)
    job='amdweight_MC'

    os.system('python {} -input {} -col {} -gamd {} -job {} {} {} -order 10 {} {} {} {} {} {} {} {} | tee -a {}.{}-reweight_MC.log'.format( 
              self.pyrw, in_file, self.col, self.gamd, job, self.binz, self.T,
              self.xdim, self.cut, self.emax, xlab, self.dpi,
              self.c_step, self.smo, self.pwd,   name, self.dimn))

    os.system('mv -v pmf-{0}.xvg pmf-{1}.{2}-MCorder10.xvg'.format(
              name+subfx, name, self.dimn))
#    os.system('xmgrace -hdevice SVG -hardcopy -nxy pmf-{0}.{1}-MCorder10.xvg -printfile pmf-{0}.{1}-MCorder10.svg'.format(name, self.dimn))


################
  def reweighting(self, rawfile):
    in_file, name, subfx, xlab = self._check_name(rawfile)
    job='amdweight'

    os.system('python {} -input {} -col {} -gamd {} -job {} {} {} {} {} {} {} {} {} {} {} | tee -a {}.{}-reweight.log'.format( 
              self.pyrw, in_file, self.col, self.gamd, job, self.binz, self.T,
              self.xdim, self.cut, self.emax, xlab, self.dpi,
              self.c_step, self.smo, self.pwd,    name, self.dimn))

    os.system('mv -v pmf-{0}.xvg pmf-{1}.{2}-reweight.xvg'.format(
              name+subfx, name, self.dimn))
    os.system('rm weights-{}.xvg'.format(name+subfx))
#    os.system('xmgrace -hdevice SVG -hardcopy -nxy pmf-{0}.{1}-reweight.xvg -printfile pmf-{0}.{1}-reweight.svg'.format(name, self.dimn))


################
  def no_weighting(self, rawfile):
    in_file, name, subfx, xlab = self._check_name(rawfile)
    job='noweight'

    os.system('python {} -input {} -col {} -gamd {} -job {} {} {} {} {} {} {} {} {} {} {} | tee -a {}.{}-no_weight.log'.format( 
              self.pyrw, in_file, self.col, self.gamd, job, self.binz, self.T,
              self.xdim, self.cut, self.emax, xlab, self.dpi,
              self.c_step, self.smo, self.pwd,   name, self.dimn))

    os.system('mv -v pmf-{0}.xvg pmf-{1}.{2}-noweight.xvg'.format(
              name+subfx, name, self.dimn))
#    os.system('xmgrace -hdevice SVG -hardcopy -nxy pmf-{0}.{1}-noweight.xvg -printfile pmf-{0}.{1}-noweight.svg'.format(name, self.dimn))


################
  def amd_time(self, rawfile):
    in_file, name, subfx, xlab = self._check_name(rawfile)
    job='amd_time'

    os.system('python {} -input {} -col {} -gamd {} -job {} {} {} | tee -a {}.{}-amd_time.log'.format( 
              self.pyrw, in_file, self.gamd, job, self.binz, self.T, self.pwd,
              name, self.dimn))


###########################################################################

def cmdparse():
  parser = ArgumentParser(description="command line arguments")

  parser.add_argument("-inp", dest="input", required=True, 
                      help="list of input files", metavar="<input list>")
  parser.add_argument("-gamd", dest="gamd", required=True,
                      help="GAMD weight file", metavar="<GAMD weight file>")

  parser.add_argument('-col', dest='col', required=False,
                      help='column number to be read (def: 2)', metavar='<column>')

  parser.add_argument('-pwd', dest='pwd', required=False,
                      help="Path of working directory, beware of characters' \()' (def: ./)", metavar='<work directory>')
  parser.add_argument('-dir', dest='dir', required=False,
                      help='path to input files only, on top of -pwd (def: ./)', metavar='<data>')
  parser.add_argument("-bin", dest="bin", required=False,
                      help="General bin size (def: 0.5)", metavar="<bin size>")
  parser.add_argument("-Emax", dest="Emax", required=False,
                      help="Max energy above 0 kcal/mol (def: 4)", metavar="<Emax>")

  parser.add_argument("-temp", dest="T", required=False, 
                      help="Temperature K (def: 310)",metavar="<temperature>")
  parser.add_argument("-Xdim", dest="Xdim", required=False, nargs="+", 
                      help="Data range in X-dimension (def: auto)", metavar="<Xmin Xmax>")

  parser.add_argument("-histcut", dest="cut", required=False,
                      help="Data cutoff (def: 25)", metavar="<histo cutoff>")
  parser.add_argument('-contour', dest='c_step', required=False,
                      help='Fig Contour per integer (def: 4)', metavar='<contour>')
  parser.add_argument('-smooth', dest='smooth', required=False,
                      help='Contour smoothening (def: None | Sug: 1.15)', metavar='<smooth>')
  parser.add_argument('-dpi', dest='dpi', required=False,
                      help='Fig Resolution (def: 200)', metavar='<dpi>')

  args = parser.parse_args()
  return(args)


##########################################################################

if __name__ == '__main__':
    main()

