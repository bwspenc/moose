#!/bin/sh
mkdir elas_ref
cd elas_ref
mpiexec -n 4 ~/gitproj/githerd/moose/modules/tensor_mechanics/tensor_mechanics-opt -i ../planestress_elasticity_weak.i Mesh/fmg/file=platehole_04_q4.e Outputs/file_base=elastic_weak_p1p1_04
mpiexec -n 4 ~/gitproj/githerd/moose/modules/tensor_mechanics/tensor_mechanics-opt -i ../planestress_elasticity_weak.i Mesh/fmg/file=platehole_02_q4.e Outputs/file_base=elastic_weak_p1p1_02
mpiexec -n 4 ~/gitproj/githerd/moose/modules/tensor_mechanics/tensor_mechanics-opt -i ../planestress_elasticity_weak.i Mesh/fmg/file=platehole_01_q4.e Outputs/file_base=elastic_weak_p1p1_01

../contour_plot_variable.py -b elastic_weak_p1p1_04 -v vonmises_stress
../contour_plot_variable.py -b elastic_weak_p1p1_02 -v vonmises_stress
../contour_plot_variable.py -b elastic_weak_p1p1_01 -v vonmises_stress
../contour_plot_variable.py -b elastic_weak_p1p1_04 -v stress_zz
../contour_plot_variable.py -b elastic_weak_p1p1_02 -v stress_zz
../contour_plot_variable.py -b elastic_weak_p1p1_01 -v stress_zz


mpiexec -n 4 ~/gitproj/githerd/moose/modules/tensor_mechanics/tensor_mechanics-opt -i ../planestress_elasticity_weak.i Mesh/fmg/file=platehole_04_q4.e Outputs/file_base=elastic_weak_p1p0_04 Variables/nl_strain_zz/family=MONOMIAL Variables/nl_strain_zz/order=CONSTANT
mpiexec -n 4 ~/gitproj/githerd/moose/modules/tensor_mechanics/tensor_mechanics-opt -i ../planestress_elasticity_weak.i Mesh/fmg/file=platehole_02_q4.e Outputs/file_base=elastic_weak_p1p0_02 Variables/nl_strain_zz/family=MONOMIAL Variables/nl_strain_zz/order=CONSTANT
mpiexec -n 4 ~/gitproj/githerd/moose/modules/tensor_mechanics/tensor_mechanics-opt -i ../planestress_elasticity_weak.i Mesh/fmg/file=platehole_01_q4.e Outputs/file_base=elastic_weak_p1p0_01 Variables/nl_strain_zz/family=MONOMIAL Variables/nl_strain_zz/order=CONSTANT

../contour_plot_variable.py -b elastic_weak_p1p0_04 -v vonmises_stress
../contour_plot_variable.py -b elastic_weak_p1p0_02 -v vonmises_stress
../contour_plot_variable.py -b elastic_weak_p1p0_01 -v vonmises_stress
../contour_plot_variable.py -b elastic_weak_p1p0_04 -v stress_zz
../contour_plot_variable.py -b elastic_weak_p1p0_02 -v stress_zz
../contour_plot_variable.py -b elastic_weak_p1p0_01 -v stress_zz


mpiexec -n 4 ~/gitproj/githerd/moose/modules/tensor_mechanics/tensor_mechanics-opt -i ../planestress_elasticity_weak.i Mesh/fmg/file=platehole_04_q9.e Outputs/file_base=elastic_weak_p2p2_04 GlobalParams/order=SECOND
mpiexec -n 4 ~/gitproj/githerd/moose/modules/tensor_mechanics/tensor_mechanics-opt -i ../planestress_elasticity_weak.i Mesh/fmg/file=platehole_02_q9.e Outputs/file_base=elastic_weak_p2p2_02 GlobalParams/order=SECOND
mpiexec -n 4 ~/gitproj/githerd/moose/modules/tensor_mechanics/tensor_mechanics-opt -i ../planestress_elasticity_weak.i Mesh/fmg/file=platehole_01_q9.e Outputs/file_base=elastic_weak_p2p2_01 GlobalParams/order=SECOND

../contour_plot_variable.py -b elastic_weak_p2p2_04 -v vonmises_stress
../contour_plot_variable.py -b elastic_weak_p2p2_02 -v vonmises_stress
../contour_plot_variable.py -b elastic_weak_p2p2_01 -v vonmises_stress
../contour_plot_variable.py -b elastic_weak_p2p2_04 -v stress_zz
../contour_plot_variable.py -b elastic_weak_p2p2_02 -v stress_zz
../contour_plot_variable.py -b elastic_weak_p2p2_01 -v stress_zz


mpiexec -n 4 ~/gitproj/githerd/moose/modules/tensor_mechanics/tensor_mechanics-opt -i ../planestress_elasticity_weak.i Mesh/fmg/file=platehole_04_q9.e Outputs/file_base=elastic_weak_p2p1_04 GlobalParams/order=SECOND Variables/nl_strain_zz/order=FIRST
mpiexec -n 4 ~/gitproj/githerd/moose/modules/tensor_mechanics/tensor_mechanics-opt -i ../planestress_elasticity_weak.i Mesh/fmg/file=platehole_02_q9.e Outputs/file_base=elastic_weak_p2p1_02 GlobalParams/order=SECOND Variables/nl_strain_zz/order=FIRST
mpiexec -n 4 ~/gitproj/githerd/moose/modules/tensor_mechanics/tensor_mechanics-opt -i ../planestress_elasticity_weak.i Mesh/fmg/file=platehole_01_q9.e Outputs/file_base=elastic_weak_p2p1_01 GlobalParams/order=SECOND Variables/nl_strain_zz/order=FIRST

../contour_plot_variable.py -b elastic_weak_p2p1_04 -v vonmises_stress
../contour_plot_variable.py -b elastic_weak_p2p1_02 -v vonmises_stress
../contour_plot_variable.py -b elastic_weak_p2p1_01 -v vonmises_stress
../contour_plot_variable.py -b elastic_weak_p2p1_04 -v stress_zz
../contour_plot_variable.py -b elastic_weak_p2p1_02 -v stress_zz
../contour_plot_variable.py -b elastic_weak_p2p1_01 -v stress_zz

mpiexec -n 4 ~/gitproj/githerd/moose/modules/tensor_mechanics/tensor_mechanics-opt -i ../planestress_elasticity_weak.i Mesh/fmg/file=platehole_04_q9.e Outputs/file_base=elastic_weak_p2p0_04 GlobalParams/order=SECOND Variables/nl_strain_zz/family=MONOMIAL Variables/nl_strain_zz/order=CONSTANT
mpiexec -n 4 ~/gitproj/githerd/moose/modules/tensor_mechanics/tensor_mechanics-opt -i ../planestress_elasticity_weak.i Mesh/fmg/file=platehole_02_q9.e Outputs/file_base=elastic_weak_p2p0_02 GlobalParams/order=SECOND Variables/nl_strain_zz/family=MONOMIAL Variables/nl_strain_zz/order=CONSTANT
mpiexec -n 4 ~/gitproj/githerd/moose/modules/tensor_mechanics/tensor_mechanics-opt -i ../planestress_elasticity_weak.i Mesh/fmg/file=platehole_01_q9.e Outputs/file_base=elastic_weak_p2p0_01 GlobalParams/order=SECOND Variables/nl_strain_zz/family=MONOMIAL Variables/nl_strain_zz/order=CONSTANT

../contour_plot_variable.py -b elastic_weak_p2p0_04 -v vonmises_stress
../contour_plot_variable.py -b elastic_weak_p2p0_02 -v vonmises_stress
../contour_plot_variable.py -b elastic_weak_p2p0_01 -v vonmises_stress
../contour_plot_variable.py -b elastic_weak_p2p0_04 -v stress_zz
../contour_plot_variable.py -b elastic_weak_p2p0_02 -v stress_zz
../contour_plot_variable.py -b elastic_weak_p2p0_01 -v stress_zz
