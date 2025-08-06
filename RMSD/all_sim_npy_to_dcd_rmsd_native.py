#!/usr/bin/env python

import os.path
import glob
import subprocess

dcds = glob.glob('/users/ssyjb1/RNA-ML/RNA_dihe/Simulation_openmm/Simulation_new/Simulation_cos/Simulation_3_30_36_34567_w=12_320/run_all/run0**/md.dcd')

for dcd in dcds:
    d, ext = os.path.splitext(dcd)
    rmsd = open(d + '.rmsd', 'w')
    subprocess.run(['/users/ssyjb1/python/lop/bestfit_dcd_by_pdb_part_onlyRMSD.py', './Simulation_new/T2HP_folded_30.pdb', dcd, '1', '30'], stdout=rmsd)
    rmsd.close()

