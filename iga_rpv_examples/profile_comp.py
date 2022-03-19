#!/usr/bin/env python

import pandas
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import sys

cases=['1d', '1d_coarse', '3d_lagrange', '3d_berenstein', '2d_lagrange']


names={'3d_berenstein': '3D 2nd order Berenstein, coarse',
        '1d': '1D 2nd order Lagrange, fine',
        '1d_coarse': '1D 2nd order Lagrange, coarse',
        '3d_lagrange': '3D 2nd order Lagrange, coarse',
        '2d_lagrange': '2D Cartesian, 2nd order Lagrange, fine'}

dataframes={}
for case in cases:
  dataframes[case]=pandas.read_csv(case+'/csv/out_hoop_stress_base_0_0_0020.csv')

f1 = plt.figure(1, figsize=(6,4.5))
ax = plt.gca()

for case in cases:
  dataframes[case].plot(ax=ax, x='id', y='hoop_stress_base', label=names[case])

ax.set_xlabel('Radial Position (m)')
ax.set_ylabel('Hoop Stress (Pa)')
plt.savefig('profile_comp.pdf', bbox_inches='tight', pad_inches=0.1)
