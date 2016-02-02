#!/usr/bin/env python

import os, re, sys

first = True
# name of the executable file
exename = "/home/uccoaku/multivarM_perturb/src/xmatfluct_two_trait_sin"

# list of values for the parameter dictionary
par_dict_values = [ "c", "mu_g", "sdmu_g", "mu_m", "sdmu_m", "phi","m11","m12","m21","m22","sigma_p","rate1","rate2","rate1ptb","rate2ptb","phiptb","int1ptb","int2ptb","diagonal_only" ]

# script that analyzes individual-based simulations 
# of multivariate maternal effects in perturbed envts

# also allows simulations that have not been finished
# and for which no parameter output has been given
# to be processed, by getting the parameters from
# the main jobfile

def analyze_parameters(lines,first=False):

    pars = {}

    for line in lines:
        mobj = re.match("(.*):;(.*);",line)
        if mobj != None:
            pars[mobj.group(1)] = mobj.group(2)

    return(pars)

# process jobfile's line of parameters
# to store it in a dict
def process_parline(line):

    global par_dict_values

    line_values = re.split("\s+", line)

    if len(line_values) != len(par_dict_values):
        raise 

    dict = {}

    for i in range(0,len(par_dict_values)):
        dict[par_dict_values[i]] = str(line_values[i])

    return dict


# find the parameters from the jobfile
def find_params_jobfile(filename):

    global exename

    # find out the number of this simulation
    # within the core folder
    dirfiles = os.listdir(os.path.dirname(filename))
    dirfiles_sim = []

    # sort out the matching files
    for dirfile in dirfiles:
        if re.match("sim_.*",dirfile) != None:
            dirfiles_sim.append(dirfile)

    dirfiles_sim.sort()

    # find out position in list:
    position = dirfiles_sim.index(os.path.basename(filename))

#    print("dirfiles: " + " ".join(dirfiles_sim))

#    print("position " + str(position))

    # now find the corresponding parameter string in the jobfile
    mo_core = re.match("(.*)core_(\d+).*",filename)
    dir_number = mo_core.group(2)
    hpcjob_dir = mo_core.group(1)

    # remove trailing slashes from hpcjob_dir
    hpcjob_dir = re.sub("\/$","",hpcjob_dir)

    # now open the corresponding hpcfile
    f = open(hpcjob_dir + "/hpcjob_" + dir_number + ".qsub")
    fl = f.readlines()
    f.close()

    job_dict = {}

    job_ctr = 0

    # loop through the lines and collect jobs
    for line in fl:

        # match a line that contains a call to an executable
        mo_job = re.match(re.compile("^" + exename + "\s(.*)$"),line.strip())

        # store the parameters
        if mo_job != None:

            if job_ctr == position:

                # process the parameters
                job_dict = process_parline(mo_job.group(1))
                break;

            job_ctr += 1

    return job_dict


def analyze_file(filename):

    global first;

    # open file; read data
    f = open(filename)
    fl = f.readlines()
    f.close

    if len(fl) < 3:
        return


    flhead = fl[0]

    lc = len(fl)
    parline = -1

    linerange = range(0,lc)
    linerange.reverse()

    # search the parameter line
    for lineno in linerange:

        if re.match("^type",fl[lineno]) != None:
            parline = lineno
            break

    if parline != 1:
        parameters = find_params_jobfile(filename)
        #parameters = analyze_parameters(fl[parline:])
        datvals = fl[parline - 3].strip()
    else:
        # no parameters; try to look it up from the jobfile
        parameters = find_params_jobfile(filename)
        datvals = fl[-1].strip()

    if first:
        print ";".join(parameters.keys()) + ";" + flhead.strip() + "file"
        first = False

    print ";".join(parameters.values()) + ";" + datvals + filename


def visit(arg, dirname, names):
    for name in names:

        if re.match("sim_.*",name) != None:

#            print dirname + "/" + name
            data = analyze_file(dirname + "/" + name)



os.path.walk(sys.argv[1], visit, None)
