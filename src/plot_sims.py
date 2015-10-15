#!/usr/bin/env python

import sys, re, os.path, csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from pylab import *

params={'axes.linewidth' : .5}
rcParams.update(params)

# get the filename from the command line
filename = sys.argv[1]

# get parameter line
#parline = -1
#f = open(filename)
#fl = f.readlines()
#f.close()
#for i in range(0,len(fl)):
#    if re.match("^type",fl[i]):
#        parline = i
#        break


# read in the csv file
dictdat = csv.DictReader(open(filename,"r"), delimiter=";")

# process the parameters at the end of the file
def process_params(dictionary, rowctr):

    fo = open(filename,"r")
    fl = fo.readlines()

    params = {};

    for line in fl[rowctr:]:
        if line.strip() != "":
            splitted = line.strip().split(";")
            params[splitted[0]] = splitted[1]

    return params



# data can now only accessed through looping row by row
# whereas we want lists of each column
# this function does that
def get_csvdict_by_column(the_raw_dict):

    # initialize a empty dict to contain the data
    by_column_dat = {}

    rowctr = 0

    # loop through the rows of the csv file and
    # put data in the dictionary
    for row in dictdat:

        if rowctr == 0:
            for key in row.keys():
                by_column_dat[key] = []

        for key, val in row.iteritems():

            if val == "type:" :
                params = process_params(dictdat,rowctr+1)

                return (params,by_column_dat)
                
            elif key != "" and val != None and val != "":
                by_column_dat[key].append(float(val))

        rowctr += 1

    

params,histdat = get_csvdict_by_column(dictdat)

# generate the figure

# initialize and specify size 
fig = plt.figure(figsize=(10,20))

num_rows = 4

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

plt.subplot(num_rows,1,3)
plt.plot(histdat["generation"],histdat["meanphen1"],'#129aff',
        histdat["generation"],histdat["meanphen2"],'#a60090',
        linewidth=1)
plt.tick_params(axis='x',which='both',bottom='on',top='on',labelbottom='off')
plt.ylabel(r'phenotype, $\bar{z}_{i}$')
plt.legend((r'$\bar{z}_{1}$',r'$\bar{z}_{2}$'))

plt.subplot(num_rows,1,4)
plt.plot(histdat["generation"],histdat["zopt1"],'b',
        histdat["generation"],histdat["zopt2"],'r',
        linewidth=1)
plt.tick_params(axis='x',which='both',bottom='on',top='on',labelbottom='off')
plt.ylabel(r'optimum, $\theta_{i}$')
plt.legend((r'$\theta_{1}$',r'$\theta_{2}$'))

graphname = os.path.dirname(filename)
if graphname != '':
    graphname += "/"
graphname += "graph_" + os.path.basename(filename) + ".pdf"

plt.savefig(graphname,format="pdf")
