#!/usr/bin/env python
import chigger
import argparse

parser = argparse.ArgumentParser(description='Plot a single specified contour variable')
parser.add_argument("-v", "--variable", type=str, help="Name of the variable to plot", required=True)
parser.add_argument("-b", "--basename", type=str, help="Base name of the .e file to plot", required=True)

args = parser.parse_args()

reader = chigger.exodus.ExodusReader(args.basename+'.e')
result = chigger.exodus.ExodusResult(reader, variable=args.variable, representation='surface')
cbar = chigger.exodus.ExodusColorBar(result, location='left', colorbar_origin = [0.35, 0.17, 0], width=0.015, length=0.2, notation='scientific', precision=2)
cbar.setOptions('primary', font_color=[0, 0, 0], num_ticks=6, font_size=16)

window = chigger.RenderWindow(result, cbar, size=[2000,2000], background=[1, 1, 1], test=True)
window.write(args.basename+'_'+args.variable+'.png')
window.start()
