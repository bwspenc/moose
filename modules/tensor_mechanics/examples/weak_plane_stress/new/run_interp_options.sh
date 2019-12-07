#!/bin/sh
mkdir p1p1
cd p1p1
mpiexec -n 4 ~/gitproj/githerd/moose/modules/tensor_mechanics/tensor_mechanics-opt -i ../planestress_elasticity_weak.i Mesh/fmg/file=platehole_02_q4.e Outputs/file_base=elastic_weak_p1p1
mpiexec -n 4 ~/gitproj/githerd/moose/modules/tensor_mechanics/tensor_mechanics-opt -i ../planestress_elasticity_3d.i   Mesh/fmg/file=platehole_02_q4.e Outputs/file_base=elastic_3d_p1p1
mpiexec -n 4 ~/gitproj/githerd/moose/modules/tensor_mechanics/tensor_mechanics-opt -i ../planestress_plasticity_weak.i Mesh/fmg/file=platehole_02_q4.e Outputs/file_base=plastic_weak_p1p1 GlobalParams/volumetric_locking_correction=true 
mpiexec -n 4 ~/gitproj/githerd/moose/modules/tensor_mechanics/tensor_mechanics-opt -i ../planestress_plasticity_3d.i   Mesh/fmg/file=platehole_02_q4.e Outputs/file_base=plastic_3d_p1p1 GlobalParams/volumetric_locking_correction=true

../contour_plot_variable.py -b elastic_weak_p1p1 -v vonmises_stress
../contour_plot_variable.py -b elastic_weak_p1p1 -v stress_zz

../contour_plot_variable.py -b elastic_3d_p1p1 -v vonmises_stress
../contour_plot_variable.py -b elastic_3d_p1p1 -v stress_zz

../contour_plot_variable.py -b plastic_weak_p1p1 -v vonmises_stress
../contour_plot_variable.py -b plastic_weak_p1p1 -v plastic_strain_eff
../contour_plot_variable.py -b plastic_weak_p1p1 -v stress_zz

../contour_plot_variable.py -b plastic_3d_p1p1 -v vonmises_stress
../contour_plot_variable.py -b plastic_3d_p1p1 -v plastic_strain_eff
../contour_plot_variable.py -b plastic_3d_p1p1 -v stress_zz

cd ..
mkdir p1p0
cd p1p0
mpiexec -n 4 ~/gitproj/githerd/moose/modules/tensor_mechanics/tensor_mechanics-opt -i ../planestress_elasticity_weak.i Mesh/fmg/file=platehole_02_q4.e Outputs/file_base=elastic_weak_p1p0 Variables/nl_strain_zz/family=MONOMIAL Variables/nl_strain_zz/order=CONSTANT
mpiexec -n 4 ~/gitproj/githerd/moose/modules/tensor_mechanics/tensor_mechanics-opt -i ../planestress_plasticity_weak.i Mesh/fmg/file=platehole_02_q4.e Outputs/file_base=plastic_weak_p1p0 GlobalParams/volumetric_locking_correction=true Variables/nl_strain_zz/family=MONOMIAL Variables/nl_strain_zz/order=CONSTANT

../contour_plot_variable.py -b elastic_weak_p1p0 -v vonmises_stress
../contour_plot_variable.py -b elastic_weak_p1p0 -v stress_zz

../contour_plot_variable.py -b plastic_weak_p1p0 -v vonmises_stress
../contour_plot_variable.py -b plastic_weak_p1p0 -v plastic_strain_eff
../contour_plot_variable.py -b plastic_weak_p1p0 -v stress_zz

cd ..
mkdir p2p2
cd p2p2
mpiexec -n 4 ~/gitproj/githerd/moose/modules/tensor_mechanics/tensor_mechanics-opt -i ../planestress_elasticity_weak.i Mesh/fmg/file=platehole_02_q9.e Outputs/file_base=elastic_weak_p2p2 GlobalParams/order=SECOND
mpiexec -n 4 ~/gitproj/githerd/moose/modules/tensor_mechanics/tensor_mechanics-opt -i ../planestress_elasticity_3d.i   Mesh/fmg/file=platehole_02_q9.e Outputs/file_base=elastic_3d_p2p2 GlobalParams/order=SECOND
mpiexec -n 4 ~/gitproj/githerd/moose/modules/tensor_mechanics/tensor_mechanics-opt -i ../planestress_plasticity_weak.i Mesh/fmg/file=platehole_02_q9.e Outputs/file_base=plastic_weak_p2p2 GlobalParams/order=SECOND
mpiexec -n 4 ~/gitproj/githerd/moose/modules/tensor_mechanics/tensor_mechanics-opt -i ../planestress_plasticity_3d.i   Mesh/fmg/file=platehole_02_q9.e Outputs/file_base=plastic_3d_p2p2 GlobalParams/order=SECOND

../contour_plot_variable.py -b elastic_weak_p2p2 -v vonmises_stress
../contour_plot_variable.py -b elastic_weak_p2p2 -v stress_zz

../contour_plot_variable.py -b elastic_3d_p2p2 -v vonmises_stress
../contour_plot_variable.py -b elastic_3d_p2p2 -v stress_zz

../contour_plot_variable.py -b plastic_weak_p2p2 -v vonmises_stress
../contour_plot_variable.py -b plastic_weak_p2p2 -v plastic_strain_eff
../contour_plot_variable.py -b plastic_weak_p2p2 -v stress_zz

../contour_plot_variable.py -b plastic_3d_p2p2 -v vonmises_stress
../contour_plot_variable.py -b plastic_3d_p2p2 -v plastic_strain_eff
../contour_plot_variable.py -b plastic_3d_p2p2 -v stress_zz

cd ..
mkdir p2p1
cd p2p1
mpiexec -n 4 ~/gitproj/githerd/moose/modules/tensor_mechanics/tensor_mechanics-opt -i ../planestress_elasticity_weak.i Mesh/fmg/file=platehole_02_q9.e Outputs/file_base=elastic_weak_p2p1 GlobalParams/order=SECOND Variables/nl_strain_zz/order=FIRST
mpiexec -n 4 ~/gitproj/githerd/moose/modules/tensor_mechanics/tensor_mechanics-opt -i ../planestress_plasticity_weak.i Mesh/fmg/file=platehole_02_q9.e Outputs/file_base=plastic_weak_p2p1 GlobalParams/order=SECOND Variables/nl_strain_zz/order=FIRST

../contour_plot_variable.py -b elastic_weak_p2p1 -v vonmises_stress
../contour_plot_variable.py -b elastic_weak_p2p1 -v stress_zz

../contour_plot_variable.py -b plastic_weak_p2p1 -v vonmises_stress
../contour_plot_variable.py -b plastic_weak_p2p1 -v plastic_strain_eff
../contour_plot_variable.py -b plastic_weak_p2p1 -v stress_zz

cd ..
mkdir p2p0
cd p2p0
mpiexec -n 4 ~/gitproj/githerd/moose/modules/tensor_mechanics/tensor_mechanics-opt -i ../planestress_elasticity_weak.i Mesh/fmg/file=platehole_02_q9.e Outputs/file_base=elastic_weak_p2p0 GlobalParams/order=SECOND Variables/nl_strain_zz/family=MONOMIAL Variables/nl_strain_zz/order=CONSTANT
mpiexec -n 4 ~/gitproj/githerd/moose/modules/tensor_mechanics/tensor_mechanics-opt -i ../planestress_plasticity_weak.i Mesh/fmg/file=platehole_02_q9.e Outputs/file_base=plastic_weak_p2p0 GlobalParams/order=SECOND Variables/nl_strain_zz/family=MONOMIAL Variables/nl_strain_zz/order=CONSTANT

../contour_plot_variable.py -b elastic_weak_p2p0 -v vonmises_stress
../contour_plot_variable.py -b elastic_weak_p2p0 -v stress_zz

../contour_plot_variable.py -b plastic_weak_p2p0 -v vonmises_stress
../contour_plot_variable.py -b plastic_weak_p2p0 -v plastic_strain_eff
../contour_plot_variable.py -b plastic_weak_p2p0 -v stress_zz
