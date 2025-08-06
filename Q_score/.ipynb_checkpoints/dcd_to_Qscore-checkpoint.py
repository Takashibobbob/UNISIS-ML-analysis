#!/usr/bin/env python

import os
import glob
import subprocess

run_dirs = sorted(glob.glob('./run_all/run*/'))

for run_dir in run_dirs:
    input_toml = os.path.join(run_dir, 'input.toml')
    
    subprocess.run(['/users/ssyjb1/sis/sis', 'input.toml'], cwd=run_dir)
    
bps = sorted(glob.glob('./run_all/run*/md.bp'))

for bp in bps:
    base, ext = os.path.splitext(bp)
    nbp_file = base + '.nbp'

    subprocess.run([
        '/users/ssyjb1/RNA-ML/Classic_simulation/sisbp_nbp.py',
        '--ss', '/users/ssyjb1/RNA-ML/RNA_dihe/Q_score/T2HP_30.db',
        bp,
        nbp_file
    ])
