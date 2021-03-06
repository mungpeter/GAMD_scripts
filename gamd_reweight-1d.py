#!/usr/bin/env python3

import sys
import csv
import glob
import math
import scipy
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
from argparse import ArgumentParser
from scipy import special
from scipy.optimize import curve_fit
## from scipy.optimize import *

import warnings
warnings.filterwarnings('ignore', category=RuntimeWarning)
warnings.filterwarnings('ignore', category=FutureWarning)

#print('''============================================================\nPyReweighting: Python scripts used to reweight accelerated \n               and scaled molecular dynamics simulations.\nAuthors: Yinglong Miao <yinglong.miao@gmail.com>\n         Bill Sinko <wsinko@gmail.com>\nLast Update: Dec 30, 2014\nNew Update:  19.03.14   by Peter Ung\n\nCitations:\n1. Sinko W, Miao Y, de Oliveira CAF, McCammon JA (2013) Population Based Reweighting of Scaled Molecular Dynamics. The Journal of Physical Chemistry B 117(42):12759-12768.\n2. Miao Y, Sinko W, Pierce L, Bucher D, Walker RC, McCammon JA (2014) Improved reweighting of accelerated molecular dynamics simulations for free energy calculation. J Chemical Theory and Computation. 10(7): 2677-2689.\n''')

###########MAIN
def main():
## Set control parameters
    plt_figs = 0

    args = cmdlineparse()       

    data = loadfiletoarray(args.input, args.col, args.pwd)
    rows = len(data)

    weights,dV = weightparse(rows, args, args.pwd)
    
    if rows != len(weights):
        sys.exit('## ERROR: files do not match in length:\n{}\t{}\n{}\t{}'.format(
                    args.input, rows, args.weight, len(weights)))

    if args.Xdim:
        binsX, discX = assignbins(args.Xdim, args)
    else:
        max_data = np.ceil(np.amax(data))
        min_data = np.floor(np.amin(data))
#        print("(max_data,min_data) = ", max_data,min_data)
        binsX, discX = assignbins([min_data, max_data], args)  ## Default bin size

    if args.Xlab:
        xlabel = args.Xlab
    else:
        xlabel = 'RC'

##  SET MAX ENERGY FOR ALL INFINITY VALUES
    if args.Emax:
        e_max = float(args.Emax)
    else :
        e_max = 10

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

##  SET flag for Gaussian fitting of deltaV
    if args.fit:
        fit = args.fit
    else :
        fit = False	# simulation temperature
##    print("gaussian fitting:", fit)

##REWEIGHTING
    if args.job == "amdweight_CE":
        hist,newedgesX,c1,c2,c3 = reweight_CE(data,hist_min,binsX,discX,dV,T,fit)
        pmf = hist2pmf(hist, hist_min, T)
        c1  = -np.multiply(1.0/beta, c1)
        c2  = -np.multiply(1.0/beta, c2)
        c3  = -np.multiply(1.0/beta, c3)

##    	pmffile = 'pmf-comp0-'+str(args.input)+'.xvg'
##	output_pmf(pmffile,pmf,binsX, xlabel)
##   	pmffile = 'pmf-comp1-'+str(args.input)+'.xvg'
##	output_pmf(pmffile,c1,binsX, xlabel)
##   	pmffile = 'pmf-comp2-'+str(args.input)+'.xvg'
##	output_pmf(pmffile,c2,binsX, xlabel)
##   	pmffile = 'pmf-comp3-'+str(args.input)+'.xvg'
##	output_pmf(pmffile,c3,binsX, xlabel)

        c12    = np.add(c1,  c2)
        c123   = np.add(c12, c3)
        pmf_c1 = np.add(pmf, c1)
        print("\033[34mpmf_min-c1 = \033[0m{}".format(np.min(pmf_c1)))
        pmf_c1 = normalize(pmf_c1, e_max)
        pmf_c2 = np.add(pmf, c12)
        print("\033[34mpmf_min-c2 = \033[0m{}".format(np.min(pmf_c2)))
        pmf_c2 = normalize(pmf_c2, e_max)
        pmf_c3 = np.add(pmf, c123)
        print("\033[34mpmf_min-c3 = \033[0m{}".format(np.min(pmf_c3)))
        pmf_c3 = normalize(pmf_c3, e_max)
    elif args.job == "amdweight_MC":
        n        = order
        MCweight = np.zeros(len(dV))
        beta_dV  = np.multiply(dV,beta)
        for x in range(0,n+1):
            MCweight = np.add(MCweight,(np.divide(np.power(beta_dV, x), 
                                float(scipy.special.factorial(x) ))))
        weights = MCweight
        hist, newedgesX = np.histogram(data, bins = binsX, weights=weights)
        hist = prephist(hist, T, e_max)
    elif args.job == "amdweight":
        hist, newedgesX = np.histogram(data, bins = binsX, weights=weights)
        hist = prephist(hist, T, e_max)
    else:
        hist, newedgesX = np.histogram(data, bins = binsX, weights=None)
        hist = prephist(hist, T, e_max)

##SAVE FREE ENERGY DATA INTO A FILE
    if args.job == "amdweight_MC" or args.job == "amdweight" or args.job == "noweight" :
        pmffile = 'pmf-'+str(args.input)+'.xvg'
        output_pmf(pmffile,hist,binsX, xlabel)
    if args.job == "amdweight_CE" :
        hist = pmf_c1
        pmffile = 'pmf-c1-'+str(args.input)+'.xvg'
        output_pmf(pmffile,hist,binsX, xlabel)

        hist = pmf_c2
        pmffile = 'pmf-c2-'+str(args.input)+'.xvg'
        output_pmf(pmffile,hist,binsX, xlabel)

        hist = pmf_c3
        pmffile = 'pmf-c3-'+str(args.input)+'.xvg'
        output_pmf(pmffile,hist,binsX, xlabel)

##SAVE WEIGHTS
    if args.job == "amdweight_MC" or args.job == "amdweight" :
        pmffile = 'weights-'+str(args.input)+'.xvg'
        ## print("len(weights) = " + str(len(weights)) + "; len(data) = " + str(len(data)) + "; len(binsX) = " + str(len(binsX)))
        output_pmf(pmffile,weights,data, xlabel)
    if args.job == "amdweight_CE" :
        hist = np.exp(c1)
        pmffile = 'weights-c1-'+str(args.input)+'.xvg'
        output_pmf(pmffile,hist,binsX, xlabel)

        hist = np.exp(c12)
        pmffile = 'weights-c2-'+str(args.input)+'.xvg'
        output_pmf(pmffile,hist,binsX, xlabel)

        hist = np.exp(c123)
        pmffile = 'weights-c3-'+str(args.input)+'.xvg'
        output_pmf(pmffile,hist,binsX, xlabel)

    if args.job == "amd_time" or args.job == "amd_dV" :
        hist, newedgesX, binf, dV_avg, dV_std, dV_anharm, dV_mat = reweight_dV(data,hist_min,binsX,discX,dV,T)

##############
    if args.job == "amd_dV" :
        pmffile = 'dV-hist-'+str(args.input) + '.xvg'
        output_dV(pmffile,dV)

        for j in range(len(dV_avg[:])):
            nf_j = int(hist[j])
            if nf_j > 0 : 
                pmffile = 'dV-hist-'+str(args.input)+'-RC'+str('%#08.2f' % binsX[j]) + '.xvg'
                output_dV(pmffile,dV_mat[j,0:nf_j])

        alpha = anharm(dV)
        print("\033[34mAnharmonicity of all dV = \033[0m{}".format(str(alpha)))

        pmffile = 'dV-anharm-'+str(args.input)+'.xvg'
        output_dV_anharm(pmffile,binsX,dV_anharm, xlabel)

        pmffile = 'dV-stat-'+str(args.input)+'.xvg'
        output_dV_stat(pmffile,binsX,dV_avg,dV_std,dV_anharm, xlabel)

        pmffile = 'dV-mat-'+str(args.input)+'.xvg'
        output_dV_mat(pmffile,binsX,hist,dV_avg,dV_std,dV_anharm,dV_mat, xlabel)

##############
    if args.job == "amd_time":
    ##CALCULATE TOTAL SIMULATIN TIME
        time = time_amd(dV,T,dV_avg,binf)
        print("\033[34mTotal simulation time: \033[0m{}".format(time))

###PLOTTING FUNCTION FOR WEIGHTS histogram
    if plt_figs :
        [hist, edges] = np.histogram(weights, bins=100)
        width = np.absolute(np.subtract(edges[0], edges[1]))
        plt.figure(1, figsize=(11,8.5))
        plt.bar(edges[:100], hist, width=width, log=True)
        plt.yscale('log')   ###if typerror is thrown delete .matplotlib/fontList.cache  file
        plt.xticks(fontsize='18')
        plt.yticks(fontsize='18')
        plt.savefig('weights.png',bbox_inches=0)
        print("FIGURE SAVED weights.png")
        plt.show()

    print( "\nEND")

###########################################################################
# Load data files with numpy, automatically skip comment lines with '#', 
# convert input to array of columns; use specified column as input 
skip = 1    # in case too much data, slice the data array by this number
def loadfiletoarray(file, col, pwd):
    fpwd = glob.glob(pwd+'/'+file)[0]
    print('\033[34m## Reading data: \033[0m'+fpwd)
    # np.loadtxt skip comment line (#)
    loaded = np.loadtxt(fpwd, usecols=[int(col)-1])[::skip]
    print( "\033[34mDATA LOADED:\033[0m\t{0} lines\n > {1}".format(len(loaded), fpwd))
    print(loaded[1:5])
    return loaded

##########################################################################
# 
def weightparse(rows, args, pwd):
    print('\033[34m## Reading GAMD data ##\033[0m')
    print(pwd)
    print(args.weight)
    if args.job == "weighthist":
        gpwd    = glob.glob(pwd+'/'+args.weight)[0]
        gamd    = np.loadtxt(gpwd)[::skip]
        print( "\033[34mDATA LOADED:\033[0m\t{0} lines\n > {1}".format(len(gamd), gpwd))
        weights = gamd[:,0]
        dV      = np.zeros(rows)
    elif args.job == "amd_time" or args.job == "amd_dV" or args.job == "amdweight" or args.job == "amdweight_MC" or args.job == "amdweight_CE" :
        gpwd    = glob.glob(pwd+'/'+args.weight)[0]
        gamd    = np.loadtxt(gpwd)[::skip]
        print( "\033[34mDATA LOADED:\033[0m\t{0} lines\n{1}".format(len(gamd), gpwd))
        weights = np.exp(gamd[:,0])
        dV      = gamd[:,2]
    elif args.job == "noweight":
        weights = np.zeros(rows)
        weights = weights + 1
        dV      = np.zeros(rows)
    else:
        print( "ERROR JOBTYPE"+ args.job+ " NOT RECOGNIZED")
        del gamd
        del weights
        del dV
    return weights, dV

##########################################################################
def assignbins(dim, args):
    minimum = float(dim[0])
    maximum = float(dim[1])
    if args.disc:
        disc = float(args.disc)
    else :
        disc = 1
    bins = np.arange(minimum,(maximum+disc),disc)
    return bins, disc

##########################################################################
def normalize(pmf, e_max):
    pmf = pmf-np.min(pmf)  ## zero value to lowest energy state
    temphist = pmf
    #set infinity free energy values to Emax
    for x in range(len(temphist[:])):
        if np.isinf(temphist[x]):
            temphist[x] = e_max
    return temphist

##########################################################################
## Convert input data population to free energy in kcal/mol
def prephist(hist, T, e_max):
    hist = np.add(hist,0.000000000000000001)  ###so that distrib
    hist = (0.0019872036*T)*np.log(hist) ####Convert popul to free energy in Kcal/mol
    print( "PMF_min = ", -np.max(hist))
    hist = np.max(hist)-hist  ## zero value to lowest energy state
    temphist = hist
    #set infinity free energy values to Emax
    for x in range(len(temphist[:])):
        if np.isinf(temphist[x]):
            temphist[x] = e_max
    return temphist

##########################################################################
def reweight_CE(data,hist_min,binsX,discX,dV,T,fit):
    hist, newedgesX = np.histogram(data, bins = binsX, weights=None)
    hist_max = np.max(hist)

    beta  = 1.0/(0.0019872036*T)
    nf    = len(data)
    nbins = len(hist)

    c1 = np.zeros(nbins)
    c2 = np.zeros(nbins)
    c3 = np.zeros(nbins)

    binf = np.zeros(nf, dtype=int) # array for storing assigned bin of each frame
    nA = np.zeros(nbins, dtype=int) # nA is equivalent to hist here
    dV_avg  = np.zeros(nbins)
    dV_avg2 = np.zeros(nbins)
    dV_avg3 = np.zeros(nbins)
    dV_std  = np.zeros(nbins)
    dV_mat  = np.zeros((nbins,hist_max)) # matrix for storing dV of each assigned 

    dV_avg_all = np.average(dV)
    dV_std_all = np.std(dV)
    print( 'dV all: avg = ', dV_avg_all, 'std = ', dV_std_all)

    diff_tol_avg = 10
    diff_tol_std = 1
    dV_binsize   = 50

    for i in range(len(data)):
        j = int((data[i]-binsX[0])/discX)
        if j >= nbins :
            j = nbins-1
        binf[i] = j
        dV_mat[j,nA[j]] = dV[i]
        nA[j] = nA[j]+1

    for j in range(nbins):
        if nA[j] >= hist_min :
            num = int(nA[j])
            atemp  = np.zeros(num)
            atemp2 = np.zeros(num)
            atemp3 = np.zeros(num)
            for k in range(num):
                atemp[k]  = dV_mat[j,k]
                atemp2[k] = dV_mat[j,k]**2
                atemp3[k] = dV_mat[j,k]**3
            if fit:
                ## calculate average/std through gaussian fitting
                hist_temp, bin_edges_temp = np.histogram(atemp, bins=dV_binsize)
                bin_centres_temp = (bin_edges_temp[:-1] + bin_edges_temp[1:])/2
                ## output original histograms
                pmffile = 'dV-hist-forFit-RC'+str('%#08.2f' % binsX[j]) + '.xvg'
                output_pmf(pmffile,hist_temp,bin_centres_temp, xlabel)
                # p0 is the initial guess for the fitting coefficients (A, mu and sigma above)
                mean = np.average(atemp)
                std  = np.std(atemp)
                p0 = [0., 1., mean, std]
                ## coeff, var_matrix = curve_fit(gauss, bin_centres_temp, hist_temp, p0=p0)
                coeff, var_matrix = curve_fit(gauss, bin_centres_temp, hist_temp, p0=p0)
                # Finally, lets get the fitting parameters, i.e. the mean and standard deviation:
                print( binsX[j], ': mean = ', coeff[2], 'sigma = ', coeff[3])
                dV_avg[j] = coeff[2]
                dV_std[j] = coeff[3]
##	              print( binsX[j], 'y0 = ', coeff[0], 'A = ', coeff[1], 'x0 = ', coeff[2], 'sigma = ', coeff[3])
##	              print( binsX[j], 'A = ', coeff[0], 'x0 = ', coeff[1], 'sigma = ', coeff[2])
                # Get the fitted curve
                hist_fit = gauss(bin_centres_temp, *coeff)
                ## output fitted histograms
                pmffile = 'dV-hist-gaussFit-RC'+str('%#08.2f' % binsX[j]) + '.xvg'
                output_pmf(pmffile,hist_fit,bin_centres_temp, xlabel)
            else:
                ## calculate average/std directly
                dV_avg[j] = np.average(atemp)
                dV_std[j] = np.std(atemp)
##                print( binsX[j], 'mean = ', dV_avg[j], 'std = ', dV_std[j])

            if np.absolute(dV_avg[j]-dV_avg_all)>diff_tol_avg or np.absolute(dV_std[j]-dV_std_all)>diff_tol_std :
                dV_avg[j] = 0
                dV_std[j] = 0
##                print( binsX[j], ': corrected mean = ', dV_avg[j], 'sigma = ', dV_std[j])

            dV_avg2[j] = np.average(atemp2)
            dV_avg3[j] = np.average(atemp3)
            del atemp
            del atemp2
            del atemp3
            c1[j] = beta*dV_avg[j]
            c2[j] = 0.5*beta**2*dV_std[j]**2
            c3[j] = (1.0/6.0)*beta**3*(dV_avg3[j]-3.0*dV_avg2[j]*dV_avg[j]+2.0*dV_avg[j]**3)
##      print( nA[j],dV_mat[j,1:10])
##	print( binsX[j],nA[j],len(atemp),dV_avg[j],dV_std[j])
    return hist, newedgesX, c1, c2, c3

##########################################################################
def reweight_dV(data,hist_min,binsX,discX,dV,T):
    hist, newedgesX = np.histogram(data, bins = binsX, weights=None)
    hist_max = np.max(hist)

    nf = len(data)
    nbins = len(hist)

    binf = np.zeros(nf, dtype=int) # array for storing assigned bin of each frame
    nA = np.zeros(nbins, dtype=int) # nA is equivalent to hist here
    dV_avg = np.zeros(nbins)
    dV_std = np.zeros(nbins)
    dV_anharm = np.zeros(nbins)
    dV_mat = np.zeros((nbins,hist_max)) # matrix for storing dV of each assigned 

    for i in range(len(data)):
        j = int((data[i]-binsX[0])/discX)
        if j >= nbins :
            j = nbins-1
        binf[i] = j
        dV_mat[j,nA[j]] = dV[i]
        nA[j] = nA[j]+1

    for j in range(nbins):
        dV_anharm[j] = 100
        if nA[j] > 0 :
            num = int(nA[j])
            atemp = np.zeros(num)
            for k in range(num):
                atemp[k] = dV_mat[j,k]
            dV_avg[j]    = np.average(atemp)
            dV_std[j]    = np.std(atemp)
            dV_anharm[j] = anharm(atemp)
            del atemp
##      print( nA[j],dV_mat[j,1:10])
##	print( binsX[j],nA[j],len(atemp),dV_avg[j],dV_std[j])
    return hist, newedgesX, binf, dV_avg, dV_std, dV_anharm, dV_mat

##########################################################################
##  Convert histogram to free energy in Kcal/mol
def hist2pmf(hist,hist_min,T):
    nbins = len(hist)
    pmf   = np.zeros(nbins)
    pmf_min = 100
    for j in range(len(hist)):
        if hist[j] >= hist_min :
            pmf[j] = -(0.0019872036*T)*np.log(hist[j])
##            print( j, pmf[j])
        if pmf_min > pmf[j] :
            pmf_min = pmf[j]
##    pmf=pmf-pmf_min  ## zero value to lowest energy state
    return pmf

##########################################################################
def output_pmf(pmffile,hist,binsX, xlabel):
    fpmf = open(pmffile, 'w')
    strpmf = '#{0} \tPMF(kcal/mol)\n\n@    xaxis  label \"{0}\"\n@    yaxis  label \"PMF(kcal/mol)\"\n@TYPE xy\n'.format(xlabel)
    fpmf.write(strpmf)
    for j in range(len(hist[:])):
        strpmf = str(binsX[j]) + ' \t' + str(hist[j]) + '\n'
        fpmf.write(strpmf)
    fpmf.closed
    return fpmf

##########################################################################
def output_dV(pmffile,dV):
    fpmf = open(pmffile, 'w')
    strpmf = '#dV \tp(dV) \n\n@    xaxis  label \"dV\"\n@    yaxis  label \"p(dV)\"\n@TYPE xy\n'
    hist_dV, bin_dV = np.histogram(dV, bins=50)
    for k in range(len(hist_dV)):
        strpmf = strpmf + str(bin_dV[k]) + ' \t' + str(hist_dV[k]) + ' \n'
        fpmf.write(strpmf)
    fpmf.closed
    return fpmf

##########################################################################
def output_dV_anharm(pmffile,binsX,dV_anharm, xlabel):
    fpmf = open(pmffile, 'w')
    strpmf = '#{0} \tdV_anharm \tError\n\n@    xaxis  label \"{0}\"\n@    yaxis  label \"dV_anmarm\"\n@TYPE xy\n'.format(xlabel)
    fpmf.write(strpmf)
    for j in range(len(dV_anharm[:])):
        strpmf = str(binsX[j]) + ' \t' + str(dV_anharm[j]) + '\n'
        fpmf.write(strpmf)
    fpmf.closed
    return fpmf

##########################################################################
def output_dV_stat(pmffile,binsX,dV_avg,dV_std,dV_anharm, xlabel):
    fpmf = open(pmffile, 'w')
    strpmf = '#{0} \tdV_avg(kcal/mol) \tError\n\n@    xaxis  label \"{0}\"\n@    yaxis  label \"dV(kcal/mol)\"\n@TYPE xydy\n'.format(xlabel)
    fpmf.write(strpmf)
    for j in range(len(dV_avg[:])):
	    strpmf = str(binsX[j]) + ' \t' + str(dV_avg[j]) + ' \t' + str(dV_std[j]) + ' \t' + str(dV_anharm[j]) + '\n'
	    fpmf.write(strpmf)
    fpmf.closed
    return fpmf

##########################################################################
def output_dV_mat(pmffile,binsX,hist,dV_avg,dV_std,dV_anharm,dV_mat, xlabel):
    fpmf = open(pmffile, 'w')
    strpmf = '#{0} \tNf \tdV_avg \tdV_std \tdV_ij \n\n@    xaxis  label \"{0}\"\n@    yaxis  label \"dV(kcal/mol)\"\n@TYPE xy\n'.format(xlabel)
    fpmf.write(strpmf)
    for j in range(len(dV_avg[:])):
        nf_j = int(hist[j])
        strpmf = str(binsX[j]) + ' \t' + str(hist[j]) + ' \t' + str(dV_avg[j]) + ' \t' + str(dV_std[j]) + ' \t' + str(dV_anharm[j])
        for k in range(int(nf_j)):
            strpmf = strpmf + ' \t' + str(dV_mat[j,k])
        strpmf = strpmf + '\n'
        fpmf.write(strpmf)
    fpmf.closed
    return fpmf

##########################################################################
def time_amd(dV,T,dV_avg,binf):
    nf    = len(dV)
    nbins = len(dV_avg)

    time = 0.0
    for i in range(nf):
        j = binf[i]
        time = time+np.exp(dV_avg[j]/(0.0019872036*T))
    return time

##########################################################################
def gauss(x, *p):
    y0, A, mu, sigma = p
    return y0+A*np.exp(-(x-mu)**2/(2.*sigma**2))

##########################################################################
def anharm(data):
#    print( "Compute anharmonicity")
    var = np.var(data)
    hist, edges = np.histogram(data, 50, density=True)
    hist = np.add(hist,0.000000000000000001)  ###so that distrib
    dx = edges[1]-edges[0]
    S1 = -1*np.trapz(np.multiply(hist, np.log(hist)),dx=dx)
    S2 = 0.5*np.log(np.add(2.00*np.pi*np.exp(1)*var,0.000000000000000001))
    alpha = S2-S1
    if np.isinf(alpha):
        alpha = 100
    return alpha

##########################################################################
def anharmND(datafull):
    print( "Performing error estimation")
    width = datafull[0,:]
    for x in range(len(width)):
        var = np.var(datafull[:,x])
        std = np.std(datafull[:,x])
        print( "variance of "+str(x)+" is : " +str(var)+" Standard Deviation:  "+str(std))
        hist, edges = np.histogram(datafull[:,x], 100, density=True)
        hist = np.add(hist,0.000000000000000001)  ###so that distrib
        dx = edges[1]-edges[0]
        S1 = -1*np.trapz(np.multiply(hist, np.log(hist)),dx=dx)
        S2 = 0.5*np.log(np.add(2.00*np.pi*np.exp(1)*var,0.000000000000000001))
        alpha = S2-S1
        print( str(x)+"dimension   S1 = "+str(S1)+"  S2 = "+str(S2)+" Alpha = "+str(alpha) ) 
    return var, std, alpha

##########################################################################
# READ datafiles and print weights
def cmdlineparse():
    parser = ArgumentParser(description="command line arguments")
    parser.add_argument("-input", dest="input", required=True, 
                        help="input file", metavar="<input file>")
    parser.add_argument('-col', dest='col', required=True,
                        help='column to be read', metavar='<column>')
    parser.add_argument("-job", dest="job", required=True, 
                        help="Reweighting method to use: <noweight>, <weighthist>, <amd_time>, <amd_dV>, <amdweight>, <amdweight_MC>, <amdweight_CE>", metavar="<Job type reweighting method>")
    parser.add_argument("-gamd", dest="weight", required=True, 
                        help="gamd file", metavar="<weight file>")
                        
    parser.add_argument("-Xdim", dest="Xdim", required=False, nargs="+", 
                        help="Data range of X-dimension (def: auto)", metavar="<Xmin Xmax>")
    parser.add_argument("-bin", dest="disc", required=False,  
                        help="Bin size", metavar="<discretization>")
    parser.add_argument("-cutoff", dest="cutoff", required=False,  
                        help="histogram cutoff (def: 10)", metavar="<cutoff>")
    parser.add_argument("-T", dest="T", required=False,  
                        help="Temperature K (def: 310)", metavar="<Temperature>")
    parser.add_argument("-Emax", dest="Emax", required=False,  
                        help="Cutoff for maximum free energy above 0 kcal/mol (def: 10)", metavar="<Emax>")
    parser.add_argument("-fit", dest="fit", required=False, 
                        help="Fit deltaV distribution (def: False)", action='store_true')
    parser.add_argument("-order", dest="order", required=False, 
                        help="Order of Maclaurin series (def: 10)", metavar="<order>")
    parser.add_argument("-Xlab", dest="Xlab", required=False,
                        help="X-axis label", metavar="<X label>")

    parser.add_argument('-dpi', dest='dpi', required=False,
                        help='Fig Resolution (def: 200)', metavar='<dpi>')
    parser.add_argument('-contour', dest='c_step', required=False,
                        help='Fig Contour per integer (def: 4)', metavar='<contour>')
    parser.add_argument('-smooth', dest='smooth', required=False,
                        help='Contour smoothening (def: None | Sug: 1.15)', metavar='<smooth>')
    parser.add_argument('-pwd', dest='pwd', required=False,
                        help='path to input data (def: ./)', metavar='<path>')

    args = parser.parse_args()
    return args

##########################################################################

if __name__ == '__main__':
    main()


##########################################################################
## modified: 19.03.14
# added separators for each definition
# updated scipy.factorial (prob v0.18 or earlier) to scipy.special.factorial
# updated nA = np.zeros(*) to np.zeros(*, dtype=int)
# updated binf = np.zeros(nf) to np.zeros(nf, dtype=int)
# updated print with print() to be python3 ready
# updated warning filter
# updated the decimal place of Boltzman const 0.001987 to 0.0019872036
# updated to print number of lines (data) input
# updated xvg labels from RC1 to xlabel name
