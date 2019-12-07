#!/bin/sh
mkdir plas_ref
cd plas_ref
mpiexec -n 4 ~/gitproj/githerd/moose/modules/tensor_mechanics/tensor_mechanics-opt -i ../planestress_plasticity_weak.i Mesh/fmg/file=platehole_04_q4.e Outputs/file_base=plastic_weak_p1p0_04 GlobalParams/volumetric_locking_correction=true Variables/nl_strain_zz/family=MONOMIAL Variables/nl_strain_zz/order=CONSTANT
mpiexec -n 4 ~/gitproj/githerd/moose/modules/tensor_mechanics/tensor_mechanics-opt -i ../planestress_plasticity_weak.i Mesh/fmg/file=platehole_02_q4.e Outputs/file_base=plastic_weak_p1p0_02 GlobalParams/volumetric_locking_correction=true Variables/nl_strain_zz/family=MONOMIAL Variables/nl_strain_zz/order=CONSTANT
mpiexec -n 4 ~/gitproj/githerd/moose/modules/tensor_mechanics/tensor_mechanics-opt -i ../planestress_plasticity_weak.i Mesh/fmg/file=platehole_01_q4.e Outputs/file_base=plastic_weak_p1p0_01 GlobalParams/volumetric_locking_correction=true Variables/nl_strain_zz/family=MONOMIAL Variables/nl_strain_zz/order=CONSTANT

../contour_plot_variable.py -b plastic_weak_p1p0_04 -v plastic_strain_eff
../contour_plot_variable.py -b plastic_weak_p1p0_02 -v plastic_strain_eff
../contour_plot_variable.py -b plastic_weak_p1p0_01 -v plastic_strain_eff
../contour_plot_variable.py -b plastic_weak_p1p0_04 -v vonmises_stress
../contour_plot_variable.py -b plastic_weak_p1p0_02 -v vonmises_stress
../contour_plot_variable.py -b plastic_weak_p1p0_01 -v vonmises_stress
../contour_plot_variable.py -b plastic_weak_p1p0_04 -v stress_zz
../contour_plot_variable.py -b plastic_weak_p1p0_02 -v stress_zz
../contour_plot_variable.py -b plastic_weak_p1p0_01 -v stress_zz
