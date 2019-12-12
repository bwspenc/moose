#!/usr/bin/env python
import chigger
import argparse
import vtk

parser = argparse.ArgumentParser(description='Plot a single specified contour variable')
parser.add_argument("-v", "--variable", type=str, help="Name of the variable to plot", required=True)
parser.add_argument("-b", "--basename", type=str, help="Base name of the .e file to plot", required=True)

args = parser.parse_args()

camera = vtk.vtkCamera()
camera.SetViewUp(0.0000, 1.0000, 0.0000)
camera.SetPosition(0.0000, 0.0000, 2.0)
camera.SetFocalPoint(0.0000, 0.0000, 0.0000)

reader = chigger.exodus.ExodusReader(args.basename+'.e')
result = chigger.exodus.ExodusResult(reader, variable=args.variable, representation='surface', camera=camera)
cbar = chigger.exodus.ExodusColorBar(result, location='right', colorbar_origin = [0.10, 0.08, 0], width=0.02, length=0.25, notation='scientific', precision=2)
cbar.setOptions('primary', font_color=[0, 0, 0], num_ticks=6, font_size=18)

window = chigger.RenderWindow(result, cbar, size=[1600,1600], background=[1, 1, 1], test=True)
window.write(args.basename+'_'+args.variable+'.png')
window.start()
