# Title: Code to generate data which performs a convergence test on the model.
# Author: Stephen Williams
# Notes: Script that can be used to produce the data needed for conv_processing_ scripts.

##

import numpy as np # Numpy for numpy
import matplotlib.pyplot as plt
from alternative_solvers.function_LW2 import LW_SPM2

# Check convergence in s-dimension.
# ds0 = 0.005
nsamples = 5 # Number of runs
dsvals = np.zeros(nsamples) # Store for ds values.
dtvals = np.zeros(nsamples) # Store for dt values.

# Set the values for ds
dsvals[0] = 0.005; dsvals[1] = 0.010
dsvals[2] = 0.020; dsvals[3] = 0.025
dsvals[4] = 0.050
dt = 0.001
## End of ds test section

filename = 'mu_test/' # Folder name for data storage.

# Execute the model for each of the set of pairs of values ds/dt.
for i,dsi in enumerate(dsvals):
    print('Entering loop ' + str(i)) # Progress update, loop start.
    N = LW_SPM2(dsi,dt,i,filename) # Model - dt run.
    print('S Loop ' + str(i) + ' Complete.') # Progress update, loop end.

##


