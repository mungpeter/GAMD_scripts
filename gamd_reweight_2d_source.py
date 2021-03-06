#!/usr/bin/env python

import sys
import csv
import math
import scipy
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy import special
from scipy import ndimage
from scipy.optimize import curve_fit
from pathos import multiprocessing
from argparse import ArgumentParser
## from scipy.optimize import *

import warnings
warnings.filterwarnings('ignore', category=RuntimeWarning)
warnings.filterwarnings('ignore', category=FutureWarning)

#print('''\n============================================================\nPyReweighting: Python scripts used to reweight accelerated\n               and scaled molecular dynamics simulations.\n\nAuthors: Yinglong Miao <yinglong.miao@gmail.com>\n         Bill Sinko <wsinko@gmail.com>\n\nLast Update: Dec 30, 2014\nNew  Update: 19.03.14	by Peter Ung\n\nCitations:\n1. Sinko W, Miao Y, de Oliveira CAF, McCammon JA (2013) Population Based Reweighting of Scaled Molecular Dynamics. The Journal of Physical Chemistry B 117(42):12759-12768.\n2. Miao Y, Sinko W, Pierce L, Bucher D, Walker RC, McCammon JA (2014) Improved reweighting of accelerated molecular dynamics simulations for free energy calculation. J Chemical Theory and Computation. 10(7): 2677-2689.\n''')

###########MAIN
def gamd_reweight_2d( args ):
## Set control parameters
    plt_figs = 1

    args = cmdlineparse()   

    if args.dir is None:
        args.dir = '.'

    data = loadfiletoarray(args.dir, args.input, args.col)

    # unique name by combining 2 different files
    ifile1 = args.input[0].split('/')[-1].split('.txt')
    ifile2 = args.input[1].split('/')[-1].split('.txt')[0].split('.')[-1]
    combiname = ifile1[0]+'.'+ifile2
    
    rows = len(data)
    weights,dV = weightparse(rows, args)

    if not args.job:
        args.job = 'amdweight_CE'

    if rows != len(weights):
        sys.exit('## ERROR: files do not match in length:\n{}\t{}\n{}\t{}'.format(
                    args.input[0], rows, args.weight, len(weights)))

    if args.Xdim:
        binsX, discX =  assignbinsX(args.Xdim, args)
    else:
        max_data = np.ceil(np.amax(data[:,0]))
        min_data = np.floor(np.amin(data[:,0]))
        binsX, discX= assignbinsX([min_data,max_data], args)  ## Default bin size
        # print( "(max_data,min_data) = ", max_data,min_data)
        # binsX, discX= assignbinsX([-180,180], args)  ## Default bin size
    if args.Ydim:
        binsY, discY = assignbinsY(args.Ydim, args)
    else:
        max_data = np.ceil(np.amax(data[:,1]))
        min_data = np.floor(np.amin(data[:,1]))
        binsY, discY = assignbinsY([min_data,max_data], args)  ## Default bin size
        # binsY, discY= assignbinsY([-180,180], args)  ## Default bin size

##  SET MAX ENERGY FOR ALL INFINITY VALUES
    if args.Emax:
        e_max = float(args.Emax)
    else :
        e_max = 5.

## SET FIGURE RESOLUTION DPI
    if args.dpi:
        dpi = int(args.dpi)
    else:
        dpi = 200

## SET FIG CONTOUR NUMBER per Emax integer
    if args.c_step:
        c_step = int(args.c_step)
    else:
        c_step = 4

## SET FIG CONTOUR SMOOTHENING LEVEL
    if args.smooth:
        smooth = float(args.smooth)
    else:
        smooth = None

##  SET HISTOGRAM CUTOFF
    if args.cutoff:
        hist_min = int(args.cutoff)
    else :
        hist_min = 10	# minimum number of configurations in one bin

##  SET ORDER of McLaurin series expansion
    if args.order:
        order = int(args.order)
    else :
        order = 10	# default

##  SET TEMPERATURE
    if args.T:
        T = float(args.T)
    else :
        T = 310	# simulation temperature
    beta = 1.0/(0.0019872036*T)

##REWEIGHTING
##  SET flag for Gaussian fitting of deltaV
    if args.fit:
        fit = args.fit
    else :
        fit = False	# simulation temperature
##    print( "gaussian fitting:", fit)

## X,Y-labels
    if args.Xlab:
        xlabel = args.Xlab
    else:
        xlabel = 'RC-1'
    if args.Ylab:
        ylabel = args.Ylab
    else:
        ylabel = 'RC-2'

##REWEIGHTING
    if args.job == "amdweight_CE":
        hist2,popul,newedgesX,newedgesY,c1,c2,c3 = reweight_CE(data,hist_min,binsX,discX,binsY,discY,dV,T,fit)
        pmf = hist2pmf2D(hist2, hist_min, T)
        c1  = -np.multiply(1.0/beta,c1)
        c2  = -np.multiply(1.0/beta,c2)
        c3  = -np.multiply(1.0/beta,c3)

        c12    = np.add(c1,c2)
        c123   = np.add(c12,c3)
        pmf_c1 = np.add(pmf, c1)
        print( "pmf_min-c1 = ", np.min(pmf_c1))
        pmf_c1 = normalize2D(pmf_c1, e_max)
        pmf_c2 = np.add(pmf, c12)
        print( "pmf_min-c2 = ", np.min(pmf_c2))
        pmf_c2 = normalize2D(pmf_c2, e_max)
        pmf_c3 = np.add(pmf, c123)
        print( "pmf_min-c3 = ", np.min(pmf_c3))
        pmf_c3 = normalize2D(pmf_c3, e_max)
    elif args.job == "amdweight_MC":
        n        = order
        MCweight = np.zeros(len(dV))
        beta_dV  = np.multiply(dV,beta)
        for x in range(0,n+1):
            MCweight = np.add(MCweight,(np.divide(np.power(beta_dV, x), 
                                np.float64(scipy.special.factorial(x)))))
        weights = MCweight
        hist2,newedgesX,newedgesY = np.histogram2d(data[:,0], data[:,1], bins = (binsX, binsY), weights=weights)
        hist2, popul = prephist(hist2, T, e_max)
    elif args.job == "amdweight":
        hist2,newedgesX,newedgesY = np.histogram2d(data[:,0], data[:,1], bins = (binsX, binsY), weights=weights)
        hist2, popul = prephist(hist2, T, e_max)
    else:
        hist2,newedgesX,newedgesY = np.histogram2d(data[:,0], data[:,1], bins = (binsX, binsY), weights=None)
        hist2, popul = prephist(hist2, T, e_max)

##SAVE FREE ENERGY DATA INTO A FILE
    if args.job == "amdweight_MC" or args.job == "amdweight" or args.job == "noweight" :
        pmffile = 'pmf-'+combiname+'.txt'+ifile1[1]+'.xvg'
        output_pmf2D(pmffile,hist2,binsX,binsY, xlabel,ylabel)
    if args.job == "amdweight_CE" :
        hist2   = pmf_c1
        pmffile = 'pmf-c1-'+combiname+'.txt'+ifile1[1]+'.xvg'
        output_pmf2D(pmffile,hist2,binsX,binsY, xlabel,ylabel)

        hist2   = pmf_c3
        pmffile = 'pmf-c3-'+combiname+'.txt'+ifile1[1]+'.xvg'
        output_pmf2D(pmffile,hist2,binsX,binsY, xlabel,ylabel)

        hist2   = pmf_c2
        pmffile = 'pmf-c2-'+combiname+'.txt'+ifile1[1]+'.xvg'
        output_pmf2D(pmffile,hist2,binsX,binsY, xlabel,ylabel)

    if args.job == "amd_dV":
        plt_figs = 0
        hist2,newedgesX,newedgesY,binfX,binfY,dV_avg,dV_std,dV_anharm,dV_mat = reweight_dV(data,hist_min,binsX,binsY,discX,discY,dV,T)

        pmffile = 'dV-hist-2D-'+combiname+'.txt'+ifile1[1]+'.xvg'
        output_dV(pmffile,dV)


        for jx in range(len(hist2[:,0])):
            for jy in range(len(hist2[0,:])):

                nf_j = int(hist2[jx,jy])
                if nf_j > 0 : 
                    pmffile = 'dV-hist-2D-'+combiname+'-RCX'+str('%#08.2f' % binsX[jx])+'-RCY'+str('%#08.2f' % binsY[jy]) + '.xvg'

        alpha = anharm(dV)
        print( "Anharmonicity of all dV = " + str(alpha))

        pmffile = 'dV-anharm-2D-'+combiname+'.txt'+ifile1[1]+'.xvg'
        output_dV_anharm2D(pmffile,binsX,binsY,dV_anharm, xlabel)

        pmffile = 'dV-stat-2D-'+combiname+'.txt'+ifile1[1]+'.xvg'
        output_dV_stat2D(pmffile,binsX,binsY,dV_avg,dV_std,dV_anharm, xlabel)

        pmffile = 'dV-mat-2D-'+combiname+'.txt'+ifile1[1]+'.xvg'
        output_dV_mat2D(pmffile,binsX,binsY,hist2,dV_avg,dV_std,dV_anharm,dV_mat, xlabel)

##########################################################################
### PLOTTING FUNCTION FOR FREE ENERGY FIGURE, smoothed with cubic spline
### interpolation
    if plt_figs :
        # plot with relative energy; *-t designation means no background 
        # color, showing only area that's been sampling in GAMD
        bar_label = 'PMF (kcal/mol)'
        PlotLandscape( hist2, weights, newedgesX, newedgesY, binsX, binsY, 
                        e_max, c_step, xlabel, ylabel, bar_label, combiname,
                        'dG',   'min',     'max',     smooth, dpi )
        PlotLandscape( hist2, weights, newedgesX, newedgesY, binsX, binsY, 
                        e_max, c_step, xlabel, ylabel, bar_label, combiname,
                        'dG-t', 'neither', 'neither', smooth, dpi )

        # plot with % population in simulation, p_max replace e_max for sidebar
        # because population can vary from very small to very large, get upper
        # limit to 5%
        if np.max(popul) < 5:
            p_max = np.rint(np.max(popul))     # round up to nearest integer
        else:
            p_max = 5
        p_step = c_step*1
        bar_label = '% Population'
        PlotLandscape(  popul, weights, newedgesX, newedgesY, binsX, binsY, 
                        p_max, p_step, xlabel, ylabel, bar_label, combiname,
                        'popul',   'min',     'max',     smooth, dpi )
        PlotLandscape(  popul, weights, newedgesX, newedgesY, binsX, binsY, 
                        p_max, p_step, xlabel, ylabel, bar_label, combiname,
                        'popul-t', 'neither', 'max', smooth, dpi )

    print("\nEND")


##########################################################################
# Load data files with numpy, automatically skip comment lines with '#', 
# convert input to array of columns; use specified column as input 
#skip = 1    # in case too much data, slice the data array by this number
def loadfiletoarray(dir, files, col, skip=1):
    
    mpi = multiprocessing.Pool(processes = multiprocessing.cpu_count()/2)
    obj = LoadData(dir=dir, col=col, skip=skip)
    Tmp = mpi.map(obj, files)
    mpi.close()
    mpi.join()

    if len(Tmp[0]) != len(Tmp[1]):
        sys.exit('## ERROR: files do not match in length:\n{0}\t{1}\n{2}\t{3}'.format(
                    files[0], len(Tmp[0]), files[1], len(Tmp[1])))

    return np.array(list(zip(Tmp[0], Tmp[1])))

##########################################################################
class LoadData(object):
    def __init__( self, dir='', col='', skip='' ):
        self.dir  = dir
        self.col  = col
        self.skip = skip

    def __call__( self, infile ):
        return self.load_txt( infile )

    # default 'comment' of np.loadtxt = '#'
    def load_txt( self, infile ):
        load = np.loadtxt(self.dir+'/'+infile, usecols=[int(self.col)-1])[::self.skip]
        print("DATA LOADED:\t{0}\t{1} lines".format(infile, len(load)))
        return load

##########################################################################
def weightparse(rows, args, skip=1):
    if args.job == "weighthist":
        data = np.loadtxt(args.weight, comments='#')[::skip]
        print( "DATA LOADED:\t{0}\t{1} lines".format(args.weight, len(data)))
        weights = data[:,0]
        dV   = np.zeros(rows)
    elif args.job == "amd_time" or args.job == "amd_dV" or args.job == "amdweight" or args.job == "amdweight_MC" or args.job == "amdweight_CE" :
        data = np.loadtxt(args.weight, comments='#')[::skip]
        print( "DATA LOADED:\t{0}\t{1} lines".format(args.weight, len(data)))
        weights = np.exp(data[:,0])
        dV = data[:,2]
    elif args.job == "noweight":
        weights = np.zeros(rows)
        weights = weights + 1
        dV      = np.zeros(rows)
    else:
        print( "ERROR JOBTYPE"+ args.job+ " NOT RECOGNIZED")
        del data
        del weights
        del dV
    return weights,dV

##########################################################################
def assignbinsX(dim, args):
    minimum = float(dim[0])
    maximum = float(dim[1])
    if args.discX:
        discX = float(args.discX)
    else :
        discX = 1
    binsX = np.arange(minimum, (maximum+discX), discX)
    return binsX, discX

##########################################################################
def assignbinsY(dim, args):
    minimum = float(dim[0])
    maximum = float(dim[1])
    if args.discY:
        discY = float(args.discY)
    else :
        discY = 1
    binsY =np.arange(minimum,(maximum+discY), discY)
    return binsY, discY

##########################################################################
def normalize2D(pmf, e_max):
    pmf = pmf-np.min(pmf)  ## zero value to lowest energy state
    temphist = pmf
    # print( "rows = ", len(temphist[:,0]))
    # print( "cols = ", len(temphist[0,:]))
    #set infinity free energy values to Emax
    for jy in range(len(temphist[0,:])):
        for jx in range(len(temphist[:,0])):
            if np.isinf(temphist[jx,jy]):
                temphist[jx,jy] = e_max
    return temphist

##########################################################################
## convert matrix of data to free energy, then to relative energy to lowest point
## numpy natural log (ln) is np.log; np.log10 is for log based 10
## R = 0.0019872036 kcal/K/mol
def prephist(hist2, T, e_max):
    popul = hist2*100/np.sum(hist2)     # normalize population to all states (% fraction)

    hist2 = np.add(hist2, 0.000000000000000001)  ###so that distrib
    hist2 = (0.0019872036*T)*np.log(hist2) ####Convert to free energy in Kcal/mol
    hist2 = np.max(hist2)-hist2  ## zero value to lowest energy state
    print(len(hist2))
    print('Max energy: {0}'.format(np.max(hist2)))
    temphist2 = hist2

    #set infinity free energy values above Emax cutoff to Emax
    for jy in range(len(temphist2[0,:])):
        for jx in range(len(temphist2[:,0])):
            if np.isinf(temphist2[jx,jy]):
                temphist2[jx,jy] = e_max
    return temphist2, popul


##########################################################################
def PlotLandscape(  hist2, weights, newedgesX, newedgesY, binsX, binsY, 
                    e_max, c_step, xlabel, ylabel, bar_label, combiname,
                    plot_type, bar_extend, plot_extend, smooth, dpi ):

    if smooth is not None:
        # using Gaussian filter based on cubic spline interpoluation to 
        # smoothen the matrix date to remove sharp edges in contour
        # x-/y-axes can be of difference scale and need independent smoothening
        Sigma  = [ (max(binsX)-min(binsX))*smooth/(len(binsX)), 
                    (max(binsY)-min(binsY))*smooth/(len(binsY)) ]
        smooth_hist2 = ndimage.filters.gaussian_filter(hist2, sigma=Sigma)
    else:
        smooth_hist2 = hist2       # no smoothening

    plt.figure(2, figsize=(11,8.5))
    colors = cm.jet

    # side bar tick, maximum = Emax value
    t_step = 1      # tick spacing per energy value
    cbar_ticks = np.linspace(0, e_max, num=(e_max*t_step)+1)
#    print('    cbar_ticks: '),
#    print(list(cbar_ticks))

    # Contour levels, set to be 'c_step' the Emax value, default is 2x
    levels = np.linspace(0, e_max, num=(e_max*c_step)+1)
#    print('plt_fig levels: '),
#    print(list(levels))

    # X- and Y-axes min and max, will be stretch to be equal
    extent = [ newedgesX[0]-1,  newedgesX[-1]+1, 
                newedgesY[-1]+1, newedgesY[0]-1   ]
#    print('plt_fig extent: '),
#    print(list(extent))

    ## create contour map; 
    # first, generate filled contour map
    plt.contourf(smooth_hist2.transpose(), origin='upper',
                    extent=extent, levels=levels, extend=plot_extend,
                    cmap=cm.get_cmap(colors, len(levels)) )

    # create colorbar instance on side based on last data input
    cbar = plt.colorbar(ticks=cbar_ticks, format=('% .1f'), 
                        extend=bar_extend, aspect=20)
    cbar.ax.set_ylabel(bar_label, rotation=270, fontsize=20, labelpad=22)

    # then impose contour lines on top
    plt.contour(smooth_hist2.transpose(), origin='upper', 
                extent=extent, levels=levels,  
                colors='black', linewidths=0.67, alpha=0.4 )

    # create gaussian heatmap to match backgound color of contour map
#    plt.imshow(smooth_hist2.transpose(), aspect='auto',
#               extent=extent, interpolation='gaussian',
#               cmap=cm.get_cmap(colors, len(levels)) )

    # Figure formatting    
    imaxes = plt.gca()
    plt.sca(cbar.ax)
    plt.clim(vmin=0, vmax=e_max)
    plt.yticks(fontsize=18)
    plt.sca(imaxes)
    axis = (min(binsX), max(binsX), min(binsY), max(binsY))
    plt.axis(axis)
    plt.xticks(size=18)
    plt.yticks(size=18)
    plt.xlabel(xlabel, fontsize=20)     # RC1
    plt.ylabel(ylabel, fontsize=20)     # RC2
    plt.savefig('2D_dG_surf.{0}.{1}.png'.format(combiname,plot_type),
                bbox_inches=0, dpi=dpi)
#    print( "FIGURE SAVED 2D_dG_surf.{0}.{1}.png".format(combiname,plot_type))

    ###PLOTTING FUNCTION FOR WEIGHTS histogram
    weight_plot = False
    if weight_plot:
        hist, edges = np.histogram(weights, bins=100)
        width = np.absolute(np.subtract(edges[0], edges[1]))
        plt.figure(1, figsize=(11,8.5))
        plt.bar(edges[:100], hist, width=width, log=True)
        plt.yscale('log')   ###if typerror is thrown delete .matplotlib/fontList.cache  file
        plt.xticks(fontsize='18')
        plt.yticks(fontsize='18')
        plt.xlabel('Frames', fontsize=20)
        plt.savefig('weights.{0}.{1}.png'.format(combiname,plot_type),
                bbox_inches=0, dpi=dpi)
#      print( "FIGURE SAVED weights.{0}.{1}.png".format(combiname,plot_type))

    plt.clf()
    plt.cla()
    plt.close('all')


##########################################################################
def reweight_CE(data,hist_min,binsX,discX,binsY,discY,dV,T,fit):

    hist1, edges1 = np.histogram(data[:,0], bins=binsX)
#    print(len(hist1), len(edges1))

    hist2, newedgesX, newedgesY = np.histogram2d(data[:,0], data[:,1], bins = (binsX, binsY), weights=None)

    hist_max = int(np.max(hist2))
##    print( np.max(hist2))

    beta   = 1.0/(0.0019872036*T)
    nf     = len(data[:,0])
    nbinsX = len(hist2[:,0])
    nbinsY = len(hist2[0,:])

    c1 = np.zeros((nbinsX,nbinsY)) 
    c2 = np.zeros((nbinsX,nbinsY)) 
    c3 = np.zeros((nbinsX,nbinsY)) 

    binfX = np.zeros(nf,dtype=int) # array for storing assigned bin of each frame
    binfY = np.zeros(nf,dtype=int) # array for storing assigned bin of each frame
    nA = np.zeros((nbinsX,nbinsY), dtype=int) # nA is equivalent to hist here
    dV_avg  = np.zeros((nbinsX,nbinsY)) 
    dV_avg2 = np.zeros((nbinsX,nbinsY)) 
    dV_avg3 = np.zeros((nbinsX,nbinsY)) 
    dV_std  = np.zeros((nbinsX,nbinsY)) 
    dV_anharm = np.zeros((nbinsX,nbinsY)) 
    dV_mat  = np.zeros((nbinsX,nbinsY,hist_max)) # matrix for storing dV of each assigned 

    dV_avg_all = np.average(dV)
    dV_std_all = np.std(dV)
    print( 'dV all: avg = ', dV_avg_all, 'std = ', dV_std_all)

    diff_tol_avg = 10
    diff_tol_std = 1
    dV_binsize   = 50

    for i in range(len(data[:,0])):
        jx = int((data[i,0]-binsX[0])/discX)
        jy = int((data[i,1]-binsY[0])/discY)
        if jx > nbinsX :
            jx = nbinsX-1
        if jy > nbinsY :
            jy = nbinsY-1
        binfX[i] = jx
        binfY[i] = jy
        # jx/jy out of bound when -Xdim/Ydim fail to cover entire data range??
        try:
            dV_mat[jx,jy,nA[jx,jy]] = dV[i]
        except IndexError:
            print(' ## dV_mat[jx,jy] out of bound: {0},{1}'.format(jx,jy))
            continue
        nA[jx,jy] = nA[jx,jy]+1

    for jx in range(nbinsX):
        for jy in range(nbinsY):
            
            dV_anharm[jx,jy] = 100

            if nA[jx,jy] >= hist_min :
                num    = int(nA[jx,jy])
                atemp  = np.zeros(num)
                atemp2 = np.zeros(num)
                atemp3 = np.zeros(num)

                for k in range(num):
                    atemp[k]  = dV_mat[jx,jy,k]
                    atemp2[k] = dV_mat[jx,jy,k]**2
                    atemp3[k] = dV_mat[jx,jy,k]**3

                dV_avg[jx,jy]    = np.average(atemp)
                dV_std[jx,jy]    = np.std(atemp)
                dV_anharm[jx,jy] = anharm(atemp)

                if np.absolute(dV_avg[jx,jy]-dV_avg_all) > diff_tol_avg or np.absolute(dV_std[jx,jy]-dV_std_all) > diff_tol_std :
                    dV_avg[jx,jy] = 0
                    dV_std[jx,jy] = 0
##	                  print( binsX[j], ': corrected mean = ', dV_avg[j], 'sigma = ', dV_std[j])

                dV_avg2[jx,jy] = np.average(atemp2)
                dV_avg3[jx,jy] = np.average(atemp3)
                del atemp
                del atemp2
                del atemp3
                c1[jx,jy] = beta*dV_avg[jx,jy]
                c2[jx,jy] = 0.5*beta**2*dV_std[jx,jy]**2
                c3[jx,jy] = (1.0/6.0)*beta**3*(dV_avg3[jx,jy]-3.0*dV_avg2[jx,jy]*dV_avg[jx,jy]+2.0*dV_avg[jx,jy]**3)

    # convert to fraction of population
    popul = hist2*100/np.sum(hist2)

    return hist2, popul, newedgesX, newedgesY, c1, c2, c3


##########################################################################
def reweight_dV(data,hist_min,binsX,binsY,discX,discY,dV,T):

    hist2, newedgesX, newedgesY = np.histogram2d(data[:,0], data[:,1], bins = (binsX, binsY), weights=None)

    hist_max = int(np.max(hist2))
##    print( np.max(hist2))

    nf     = len(data[:,0])
    nbinsX = len(hist2[:,0])
    nbinsY = len(hist2[0,:])

    binfX = np.zeros(nf,dtype=int) # array for storing assigned bin of each frame
    binfY = np.zeros(nf,dtype=int) # array for storing assigned bin of each frame
    nA = np.zeros((nbinsX,nbinsY), dtype=int) # nA is equivalent to hist here
    dV_avg = np.zeros((nbinsX,nbinsY)) 
    dV_std = np.zeros((nbinsX,nbinsY)) 
    dV_anharm = np.zeros((nbinsX,nbinsY)) 
    dV_mat = np.zeros((nbinsX,nbinsY,hist_max)) # matrix for storing dV of each assigned 

    for i in range(len(data[:,0])):
        jx = int((data[i,0]-binsX[0])/discX)
        jy = int((data[i,1]-binsY[0])/discY)
        if jx > nbinsX :
            jx = nbinsX-1
        if jy > nbinsY :
            jy = nbinsY-1
        binfX[i] = jx
        binfY[i] = jy
        # jx/jy out of bound when -Xdim/Ydim fail to cover entire data range
        try:
            dV_mat[jx,jy,nA[jx,jy]] = dV[i]
        except IndexError:
            print(' ## dV_mat[jx,jy] out of bound: {0},{1}'.format(jx,jy))
            continue
        nA[jx,jy] = nA[jx,jy]+1

    for jx in range(nbinsX):
        for jy in range(nbinsY):

            dV_anharm[jx,jy] = 100
            if nA[jx,jy] >= hist_min :
                num   = int(nA[jx,jy])
                atemp = np.zeros(num)

                for k in range(num):
                    atemp[k] = dV_mat[jx,jy,k]

                dV_avg[jx,jy]    = np.average(atemp)
                dV_std[jx,jy]    = np.std(atemp)
                dV_anharm[jx,jy] = anharm(atemp)
                del atemp
    return hist2,newedgesX,newedgesY,binfX,binfY,dV_avg,dV_std,dV_anharm,dV_mat

##########################################################################
##  Convert histogram to free energy in Kcal/mol
def hist2pmf2D(hist,hist_min,T):
    nbinsX = len(hist[:,0])
    nbinsY = len(hist[0,:])
    pmf = np.zeros((nbinsX,nbinsY))
    pmf_min = 100
    for jx in range(len(hist[:,0])):
        for jy in range(len(hist[0,:])):
            if hist[jx,jy] >= hist_min :
                pmf[jx,jy] =  -(0.0019872036*T)*np.log(hist[jx,jy])
            if pmf_min > pmf[jx,jy] :
                pmf_min = pmf[jx,jy]
##    pmf=pmf-pmf_min  ## zero value to lowest energy state
    return pmf

##########################################################################
def output_pmf2D(pmffile,hist,binsX,binsY, xlabel,ylabel):
    fpmf   = open(pmffile, 'w')
    strpmf = '#{0}\t{1}\tPMF(kcal/mol)\n\n@    xaxis  label \"{0}\"\n@    yaxis  label \"{1}\"\n@TYPE xy\n'.format(xlabel,ylabel)
    fpmf.write(strpmf)
    for jx in range(len(hist[:,0])):
        for jy in range(len(hist[0,:])):
            strpmf = str(binsX[jx]) + ' \t' + str(binsY[jy]) + ' \t' + str(hist[jx,jy]) + '\n'
            fpmf.write(strpmf)
    fpmf.closed
    return fpmf

##########################################################################
def output_dV(pmffile,dV):
    fpmf   = open(pmffile, 'w')
    strpmf = '#dV \tp(dV) \n\n@    xaxis  label \"dV\"\n@    yaxis  label \"p(dV)\"\n@TYPE xy\n'
    hist_dV, bin_dV = np.histogram(dV, bins=50)
    for k in range(len(hist_dV)):
        strpmf = strpmf + str(bin_dV[k]) + ' \t' + str(hist_dV[k]) + ' \n'
    fpmf.write(strpmf)
    fpmf.closed
    return fpmf

##########################################################################
def output_dV_anharm2D(pmffile,binsX,binsY,dV_anharm, xlabel):
    fpmf   = open(pmffile, 'w')
    strpmf = '#{0} \tdV_anharm \tError\n\n@    xaxis  label \"{0}\"\n@    yaxis  label \"dV_anmarm\"\n@TYPE xy\n'.format(xlabel)
    fpmf.write(strpmf)
    for jx in range(len(dV_anharm[:,0])):
        for jy in range(len(dV_anharm[0,:])):
            strpmf = str(binsX[jx]) + ' \t' + str(binsY[jy]) + ' \t' + str(dV_anharm[jx,jy]) + '\n'
            fpmf.write(strpmf)
    fpmf.closed
    return fpmf

##########################################################################
def output_dV_stat2D(pmffile,binsX,binsY,dV_avg,dV_std,dV_anharm, xlabel):
    fpmf   = open(pmffile, 'w')
    strpmf = '#{0} \tdV_avg(kcal/mol) \tError\n\n@    xaxis  label \"{0}\"\n@    yaxis  label \"dV(kcal/mol)\"\n@TYPE xydy\n'.format(xlabel)
    fpmf.write(strpmf)
    for jx in range(len(dV_anharm[:,0])):
        for jy in range(len(dV_anharm[0,:])):
            strpmf = str(binsX[jx]) + ' \t' + str(binsY[jy]) + ' \t' + str(dV_avg[jx,jy]) + ' \t' + str(dV_std[jx,jy]) + ' \t' + str(dV_anharm[jx,jy]) + '\n'
            fpmf.write(strpmf)
    fpmf.closed
    return fpmf

##########################################################################
def output_dV_mat2D(pmffile,binsX,binsY,hist,dV_avg,dV_std,dV_anharm,dV_mat,
                    xlabel):
    fpmf   = open(pmffile, 'w')
    strpmf = '#{0} \tNf \tdV_avg \tdV_std \tdV_ij \n\n@    xaxis  label \"{0}\"\n@    yaxis  label \"dV(kcal/mol)\"\n@TYPE xy\n'.fromat(xlabel)
    fpmf.write(strpmf)
    for jx in range(len(hist[:,0])):
        for jy in range(len(hist[0,:])):
            nf_j = int(hist[jx,jy])
            strpmf = str(binsX[jx]) + ' \t' + str(binsY[jy]) + ' \t' + str(hist[jx,jy]) + ' \t' + str(dV_avg[jx,jy]) + ' \t' + str(dV_std[jx,jy]) + ' \t' + str(dV_anharm[jx,jy])
            for k in range(int(nf_j)):
                strpmf = strpmf + ' \t' + str(dV_mat[jx,jy,k])
            strpmf = strpmf + '\n'
            fpmf.write(strpmf)
    fpmf.closed
    return fpmf

##########################################################################
def anharm(data):
#    print( "Compute anharmonicity")
    var = np.var(data)
    hist, edges = np.histogram(data, 50, density=True)
    hist = np.add(hist,0.000000000000000001)  ###so that distrib
    dx = edges[1]-edges[0]
    S1 = -1*np.trapz(np.multiply(hist, np.log(hist)),dx=dx)
    S2 = 0.5*np.log(2.00*np.pi*np.exp(1.0)*var+0.000000000000000001)
    alpha = S2-S1
    if np.isinf(alpha):
        alpha = 100
    return alpha

########################################################################
#READ datafiles and print weights
def cmdlineparse():
    parser = ArgumentParser(description="command line arguments")
    parser.add_argument("-input", dest="input", required=True, nargs="+",
                        help="Input files (Default: 2; Extra as background contour)", metavar="<X-input Y-input *Extra*>")
    parser.add_argument('-col', dest='col', required=True,
                        help='Column to be read (Usually: 2)', metavar='<column>')
    parser.add_argument("-weight", dest="weight", required=True,
                        help="GAMD weight file", metavar="<GAMD weight file>")

    parser.add_argument("-job", dest="job", required=False,
                        help="Reweighting method to use: <noweight>, <weighthist>, <amd_time>, <amd_dV>, <amdweight>, <amdweight_MC>, <amdweight_CE>: (Default: amdweight_CE)", metavar="<Job type reweighting method>")
    parser.add_argument('-dir', dest='dir', required=False,
                        help='path to input data (def: ./)', metavar='<path>')

    parser.add_argument("-Xdim", dest="Xdim", required=False, nargs="+",
                        help="Xdimensions", metavar="<Xmin Xmax >")
    parser.add_argument("-Ydim", dest="Ydim", required=False, nargs="+",
                        help="Ydimension", metavar="<Ymin Ymax >")
    parser.add_argument("-discX", dest="discX", required=False,
                        help="Discretization (bin) size in X-axis (Default: 0.5)", metavar="<discretization-X>")
    parser.add_argument("-discY", dest="discY", required=False,
                        help="Discretization (bin) size in Y-axis (Default: 0.5)", metavar="<discretization-Y>")
    parser.add_argument("-cutoff", dest="cutoff", required=False,
                        help="Histogram cutoff (Default: 10)", metavar="<cutoff>")
    parser.add_argument("-T", dest="T", required=False,
                        help="Temperature (Default: 310 K)", metavar="<Temperature>")
    parser.add_argument("-Emax", dest="Emax", required=False,
                        help="Maximum free energy (Default: 5 kcal/mol)", metavar="<Emax>")
    parser.add_argument("-fit", dest="fit", required=False,
                        help="Fit deltaV distribution", metavar="<fit>")
    parser.add_argument("-order", dest="order", required=False,
                        help="Order of Maclaurin series (Default: 10)", metavar="<order>")

    parser.add_argument("-Xlab", dest="Xlab", required=False,
                        help="2D Fig X-label (Default: RC-1)", metavar='<Xlab>')
    parser.add_argument("-Ylab", dest="Ylab", required=False,
                        help="2D Fig Y-label (Default: RC-2)", metavar='<Ylab>')
    parser.add_argument("-contour", dest="c_step", required=False,
                        help="Fig Contour step between integers (Default: 4)", metavar='<contour>')
    parser.add_argument("-smooth", dest="smooth", required=False,
                        help="Fig Contour smoothening (Def: 0| Sug: 1.15)", metavar='<smooth>')

    parser.add_argument("-dpi",  dest='dpi', required=False,
                        help='Fig Resolution (Default: 200)', metavar='<dpi>')
    args=parser.parse_args()
    return args

##########################################################################

if __name__ == '__main__':
    gamd_reweight_2d()


##########################################################################
## modified: 19.03.14
# added separators for each definition
# updated scipy.factorial (prob v0.18 or earlier) to scipy.special.factorial
# updated nA = np.zeros(*) to np.zeros(*, dtype=int)
# updated binf = np.zeros(nf) to np.zeros(nf, dtype=int)
# updated print to print() to be python3 ready
# updated warning filter
# updated np.histogram(normed=True) to density=True
# updated matplotlib pyplot.axes(ax) to pyplot.sca(ax)
# updated output png '2D_Free_Eneregy_Surface.png' to '2D_dG_surf.{}.png'
# updated Xlab and Ylab to dG_surf map png
# updated the decimal place of Boltzman const 0.001987 to 0.0019872036
# updated -input to take in 2 single-col files to use as 1 double-col file
# updated to print number of lines (data) input
# updated xvg labels from RC1/2 to xylabel names
# updated jx/jy with try to catch out of bound data if -Xdim/Ydim is too small
# updated output png resolution, format, contour etc
# updated all indentation, removed tab, to compile with python3 standard
# updated to include MD % population (fraction) figure in addition to energy
# updated data readin to use multiprocessing
# updated add data contour smoothing option
