#!/usr/bin/env python
import numpy as np
from scipy import linalg
import sys

if len(sys.argv) !=3:
    print('SCRIPT [input_coor][output]')
    sys.exit(2)

input_file = sys.argv[1]
output_file = sys.argv[2]


# Load the coordinate file
coordinates_align = np.load(input_file)

# Subtract mean to get mean-free values
mean_value = np.mean(coordinates_align,axis=0)
mean_free = coordinates_align - mean_value

# Define the time-lagged covariance function
def time_lagged_covariance(x, tau):
    T = len(x)
    num_bead = len(x[0])
    d = len(x[0][0])
    total_dim = num_bead * d
    c = np.zeros((total_dim, total_dim))

    for i in range(total_dim):
        bead_i = i // d  
        coord_i = i % d  

        for j in range(total_dim):
            bead_j = j // d  
            coord_j = j % d  

            cij = 0.
            for t in range(T - tau):
                cij += x[t, bead_i, coord_i] * x[t + tau, bead_j, coord_j]
            cij /= (T - tau - 1)  

            c[i, j] = cij  

    return c

# Calculate time-lagged covariance matrices
c = time_lagged_covariance(mean_free, 1000)
c_sym = 0.5 * (c + c.T)

c0 = time_lagged_covariance(mean_free, 0)

# Perform TICA using eigenvalue decomposition
eigvals, eigvecs = linalg.eig(c_sym, c0)

# Only keep the real part if the imaginary part is negligible
eigvals = np.real(eigvals)
eigvecs = np.real(eigvecs)

# Sort eigenvalues and eigenvectors by eigenvalues in descending order
idx = eigvals.argsort()[::-1]   
Lambda = eigvals[idx]
U = eigvecs[:, idx]

# Get the first two TICA components and reshape to (30,3) to be alligned with mean_free shape
TICA1 = U[:, 0].reshape(27,3)
TICA2 = U[:, 1].reshape(27,3)

T = len(mean_free)
#X_pc1 = np.zeros((T,))
X_tica1 = np.zeros((T,))
X_tica2 = np.zeros((T,))

# Project onto principal components
for t in range(T):
    #X_pc1[t] = mean_free[t] @ TICA1
    X_tica1[t] = np.real(mean_free[t].flatten() @ TICA1.flatten())
    X_tica2[t] = np.real(mean_free[t].flatten() @ TICA2.flatten())

# Compute histograms for PCA and TICA
tica_data = np.column_stack((TICA1, TICA2))
np.savetxt(f'{output_file}_TICA.txt', tica_data, header='TICA1 TICA2', fmt='%f')

# Compute histograms for TICA1
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
