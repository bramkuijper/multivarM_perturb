#!/usr/bin/env python 

import numpy
import math

exe = "./xmatfluct_two_trait_sin"


c = math.sqrt(0.7)
mu_g = 0.01
sdmu_g = 0.02
mu_m = [ 0.01, 0 ]
sdmu_m = 0.02
phi = [ 1 ]
m = " 0 0 0 0 "
var_p  = 0.1

# the rate of change in each environment
rate1 = list(numpy.arange(0,2*math.pi+2*math.pi/25,2*math.pi/25))
rate2 = list(numpy.arange(0,2*math.pi+2*math.pi/25,2*math.pi/25))

# vary the size of change in the mean environment 
int_size1 = [ 0, 1, 5 ]
int_size2 = [ 0, 1, 5 ]

# vary the direction of change in the mean environment (pos or neg)
direction1 = [ -1, 1 ]
direction2 = [ -1, 1 ]

ctr = 0;

for int_size_i in int_size1:
    for int_size_j in int_size2:
        for direction_i in direction1:
            for direction_j in direction2:

                # calculate the new intercepts
                int1ptb = int_size_i * direction_i
                int2ptb = int_size_j * direction_j

                for mu_mi in mu_m:
                    for rate1_i in rate1:
                        for rate2_i in rate2:
                            if rate1_i > rate2_i:
                                continue

                            for phi_i in phi:
                                phiptb = phi_i

                                ctr += 1
                                print("echo " + str(ctr))
                                print(exe + " " + str(c) + " "  + str(mu_g) + " " + str(sdmu_g) + " " + str(mu_mi) + " " + str(sdmu_m) + " " + str(phi_i) + " " \
                                        + m + " " + str(var_p) + " " + str(rate1_i) + " " + str(rate2_i) + " " + str(rate1_i) + " " + str(rate2_i) + " " \
                                        + str(phiptb) + " " + str(int1ptb) + " " + str(int2ptb) + " 0")
