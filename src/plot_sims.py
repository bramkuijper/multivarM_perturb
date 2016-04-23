#!/usr/bin/env python3

import sys, re, os.path
import pandas as pd
import matplotlib

import numpy as np
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from pylab import *

params={'axes.linewidth' : .5}
rcParams.update(params)

# get the filename from the command line
filename = sys.argv[1]

filename = sys.argv[1]

f = open(filename);
fl = f.readlines();
f.close()


parline = -1

for idx, line in enumerate(fl):
    if re.match("^type.*",line) != None:
        parline = idx - 1;
        break;

# read in the csv file
if parline > 0:
    histdat = pd.read_csv(filename, nrows=parline-3, sep=";")
else:
    histdat = pd.read_csv(filename, sep=";")

# only take every tenth generation, otherwise too much data....
histdat = histdat[histdat["generation"] % 10 == 0]


# generate the figure

# initialize and specify size 
fig = plt.figure(figsize=(10,20))

num_rows = 6


def calc_evs(row):
    tr = row["meanm11"] + row["meanm22"]
    det = row["meanm11"] * row["meanm22"] - row["meanm12"] * row["meanm21"]


    ev1 = np.nan
    ev2 = np.nan

    A = .5 * tr
    B = .5 * sqrt(abs(tr**2 - 4 * det))
    
    if tr**2 - 4 * det >= 0:
        ev1 = .5 * (tr + math.sqrt(tr**2 - 4 * det))
        ev2 = .5 * (tr - math.sqrt(tr**2 - 4 * det))
    return(pd.Series([A,B,ev1,ev2], index=['real','im','ev1','ev2']))

histdat[['real','im','ev1','ev2']] = histdat.apply(calc_evs, axis=1)

# add first subplot
plt.subplot(num_rows,1,1)
plt.plot(histdat["generation"],histdat["meang1"],'b',
        histdat["generation"],histdat["meang2"],'r',linewidth=1)
plt.tick_params(axis='x',which='both',bottom='on',top='on',labelbottom='off')
plt.ylabel(r'elevation, $\bar{a}$')
plt.legend((r'$\bar{a}_{1}$',r'$\bar{a}_{2}$'))

# add second subplot
plt.subplot(num_rows,1,2)
plt.plot(histdat["generation"],histdat["meanm11"],'c',
        histdat["generation"],histdat["meanm12"],'m',
        histdat["generation"],histdat["meanm21"],'y',
        histdat["generation"],histdat["meanm22"],'g',
        linewidth=1)
plt.tick_params(axis='x',which='both',bottom='on',top='on',labelbottom='off')
plt.ylabel(r'maternal effect, $\bar{m}_{ij}$')
plt.legend((r'$\bar{m}_{11}$',r'$\bar{m}_{12}$',r'$\bar{m}_{21}$',r'$\bar{m}_{22}$'))
plt.ylim(-1.5,1.5)

# add third subplot
plt.subplot(num_rows,1,3)
plt.plot(histdat["generation"],histdat["meanb11"],'r',
        histdat["generation"],histdat["meanb12"],'b',
        histdat["generation"],histdat["meanb21"],'m',
        histdat["generation"],histdat["meanb22"],'k',
        linewidth=1)
plt.tick_params(axis='x',which='both',bottom='on',top='on',labelbottom='off')
plt.ylabel(r'plasticity, $\bar{m}_{ij}$')
plt.legend((r'$\bar{b}_{11}$',r'$\bar{b}_{12}$',r'$\bar{b}_{21}$',r'$\bar{b}_{22}$'))
plt.ylim(-1.5,1.5)

# mean phenotypes
plt.subplot(num_rows,1,4)
plt.plot(histdat["generation"],histdat["meanphen1"],'#129aff',
        histdat["generation"],histdat["meanphen2"],'#a60090',
        linewidth=1)
plt.tick_params(axis='x',which='both',bottom='on',top='on',labelbottom='off')
plt.ylabel(r'phenotype, $\bar{z}_{i}$')
plt.legend((r'$\bar{z}_{1}$',r'$\bar{z}_{2}$'))

# optima
plt.subplot(num_rows,1,5)
plt.plot(histdat["generation"],histdat["zopt1"],'b',
        histdat["generation"],histdat["zopt2"],'r',
        linewidth=1)
plt.tick_params(axis='x',which='both',bottom='on',top='on',labelbottom='off')
plt.ylabel(r'optimum, $\theta_{i}$')
plt.legend((r'$\theta_{1}$',r'$\theta_{2}$'))


# eigenvalues M
plt.subplot(num_rows,1,6)
plt.plot(histdat["generation"],histdat["ev1"],'b',
        histdat["generation"],histdat["ev2"],'r',
        histdat["generation"],histdat["real"],'#009138',
        histdat["generation"],histdat["im"],'#9c0094',
        linewidth=1)
plt.tick_params(axis='x',which='both',bottom='on',top='on',labelbottom='off')
plt.ylabel(r'optimum, $\theta_{i}$')
plt.legend((r'$\lambda_{1}$',r'$\lambda_{2}$',r'$\mathrm{Re}$',r'$\mathrm{Im}$'))

graphname = os.path.dirname(filename)
if graphname != '':
    graphname += "/"
graphname += "graph_" + os.path.basename(filename) + ".pdf"

plt.savefig(graphname,format="pdf")
