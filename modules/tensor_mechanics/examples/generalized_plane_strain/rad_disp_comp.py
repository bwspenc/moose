#!/usr/bin/env python
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd

gps=pd.read_csv('fuel_gps_cross_section_out_disp_0001.csv')
axisym_2d=pd.read_csv('fuel_2d_axisym_out_disp_0001.csv')

ax=gps.plot(x='x', y='disp_x', label='gps')
axisym_2d.plot(ax=ax,x='x', y='disp_x', label='axisym_2d')
plt.savefig('rad_disp_gps_axisym.pdf')
