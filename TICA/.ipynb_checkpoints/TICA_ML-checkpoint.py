#!/usr/bin/env python
import numpy as np
from scipy import linalg
import sys

if len(sys.argv) !=5:
    print('SCRIPT [input_coor][mean_ref][TICA_cov][output]')
    sys.exit(2)

input_file = sys.argv[1]
mean_ref = sys.argv[2]
TICA_cov = sys.argv[3]
output_file = sys.argv[4]

# Load the coordinate file
coordinates_align = np.load(input_file)

# Subtract mean to get mean-free values
mean_value = np.loadtxt(mean_ref)
mean_free = coordinates_align - mean_value

T = len(mean_free)
#X_pc1 = np.zeros((T,))
X_tica1 = np.zeros((T,))
X_tica2 = np.zeros((T,))

# extract the TICA1 from TICA.txt
TICA1 = []
TICA2 = []
with open (TICA_cov) as f:
    for line in f:
        if line.startswith('#'):
            continue
        TICA1.append([float(line.split()[0]), float(line.split()[1]), float(line.split()[2])])
        TICA2.append([float(line.split()[3]), float(line.split()[4]), float(line.split()[5])])

TICA1 = np.array(TICA1).flatten()
TICA2 = np.array(TICA2).flatten()

# Project onto principal components
for t in range(T):
    #X_pc1[t] = mean_free[t] @ TICA1
    X_tica1[t] = np.real(mean_free[t].flatten() @ TICA1.flatten())
    X_tica2[t] = np.real(mean_free[t].flatten() @ TICA2.flatten())

# Compute histograms for TICA1
#H_tica, B_tica = np.histogram(4 * X_tica1, 50)
H_tica, B_tica = np.histogram(X_tica1, 50)
H_tica = H_tica.astype(float) / T

# compute 2D histogram
H, xedges, yedges = np.histogram2d(X_tica1, X_tica2, bins=50)
H = H.astype(float) / T

# Save the histogram and bin edges in separate files
np.savetxt(f'{output_file}_H_tica.txt', H_tica, header='H_tica', fmt='%f')
np.savetxt(f'{output_file}_B_tica.txt', B_tica, header='B_tica', fmt='%f')
np.savetxt(f'{output_file}_H.txt', H, header='H', fmt='%f')
np.savetxt(f'{output_file}_xedges.txt', xedges, header='xedges', fmt='%f')
np.savetxt(f'{output_file}_yedges.txt', yedges, header='yedges', fmt='%f')
np.savetxt(f'{output_file}_X_tica1.txt', X_tica1, header='X_tica1', fmt='%f')
np.savetxt(f'{output_file}_X_tica2.txt', X_tica2, header='X_tica2', fmt='%f')