# Test for displacement of pinched cylinder
# Ref: Figure 10 and Table 6 from Dvorkin and Bathe, Eng. Comput., Vol. 1, 1984.

# A cylinder of radius 1 m and length 2 m (along Z axis) with clamped ends
# (at z = 0 and 2 m) is pinched at mid-length by placing point loads of 10 N
# at (1, 0, 1) and (-1, 0, 1). Due to the symmetry of the problem, only 1/8th
# of the cylinder needs to be modeled.

# The normalized series solution for the displacement at the loading point is
# w = Wc E t / P = 164.24; where Wc is the displacement in m, E is the Young's
# modulus, t is the thickness and P is the point load.

# For this problem, E = 1e6 Pa, L = 2 m, R = 1 m, t = 0.01 m, P = 10 N and
# Poisson's ratio = 0.3. FEM results from different mesh discretizations are
# presented below. Only the 10x10 mesh is included as a test.

# Mesh of 1/8 cylinder |  FEM/analytical (Moose) | FEM/analytical (Dvorkin)
#                      |ratio of normalized disp.| ratio of normalized disp.
#----------------------|-------------------------|-------------------------
#     10 x 10          |          0.806          |        0.83
#     20 x 20          |          1.06           |        0.96
#     40 x 40          |          0.95           |         -
#     80 x 160         |          0.96           |         -

# The results from FEM analysis matches well with the series solution and with
# the solution presented by Dvorkin and Bathe (1984).

[Mesh]
  [./mesh]
    type = FileMeshGenerator
    file = cylinder_axischanged.e
  [../]
[]


[Variables]
  [./disp_x]
    order = FIRST
    family = LAGRANGE
  [../]
  [./disp_y]
    order = FIRST
    family = LAGRANGE
  [../]
  [./disp_z]
    order = FIRST
    family = LAGRANGE
  [../]
  [./rot_x]
    order = FIRST
    family = LAGRANGE
  [../]
  [./rot_y]
    order = FIRST
    family = LAGRANGE
  [../]
  [./rot_z]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[BCs]
  [./simply_support_x]
    type = DirichletBC
    variable = disp_x
    boundary = 'CD AB'
    value = 0.0
  [../]
  [./simply_support_y]
    type = DirichletBC
    variable = disp_y
    boundary = 'AB BC'
    value = 0.0
  [../]
  [./simply_support_z]
    type = DirichletBC
    variable = disp_z
    boundary = 'AD AB'
    value = 0.0
  [../]
  [./simply_support_rot_x]
    type = DirichletBC
    variable = rot_x
    boundary = 'AD BC AB'
    value = 0.0
  [../]
  [./simply_support_rot_y]
    type = DirichletBC
    variable = rot_y
    boundary = 'CD AD AB'
    value = 0.0
  [../]
  [./simply_support_rot_z]
    type = DirichletBC
    variable = rot_z
    boundary = 'CD AB BC'
    value = 0.0
  [../]
[]

[DiracKernels]
  # [./point1]
  #   type = ConstantPointSource
  #   variable = disp_z
  #   point = '0 0 1'
  #   value = -2.5 # P = 10
  # [../]
  [./point1]
    type = ConstantPointSource
    variable = disp_y
    point = '0 1 0'
    value = -2.5 # P = 10
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  line_search = 'none'
  petsc_options_iname = '-pc_type -sub_pc_type -sub_pc_factor_shift_type'
  petsc_options_value = 'asm      ilu          nonzero'
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-8
  dt = 1.0
  dtmin = 1.0
  end_time = 1.0
[]

[Kernels]
  [./solid_disp_x]
    type = ADStressDivergenceShell2
    block = '100'
    component = 0
    variable = disp_x
    through_thickness_order = SECOND
  [../]
  [./solid_disp_y]
    type = ADStressDivergenceShell2
    block = '100'
    component = 1
    variable = disp_y
    through_thickness_order = SECOND
  [../]
  [./solid_disp_z]
    type = ADStressDivergenceShell2
    block = '100'
    component = 2
    variable = disp_z
    through_thickness_order = SECOND
  [../]
  [./solid_rot_x]
    type = ADStressDivergenceShell2
    block = '100'
    component = 3
    variable = rot_x
    through_thickness_order = SECOND
    penalty = 1e5
  [../]
  [./solid_rot_y]
    type = ADStressDivergenceShell2
    block = '100'
    component = 4
    variable = rot_y
    through_thickness_order = SECOND
    penalty = 1e5
  [../]
  [./solid_rot_z]
    type = ADStressDivergenceShell2
    block = '100'
    component = 5
    variable = rot_z
    through_thickness_order = SECOND
    penalty = 1e5
  [../]
[]

[Materials]
  [./elasticity]
    type = ADComputeIsotropicElasticityTensorShell
    youngs_modulus = 1e6
    poissons_ratio = 0.3
    block = '100'
    through_thickness_order = SECOND
  [../]
  [./strain]
    type = ADComputeIncrementalShellStrain2
    block = '100'
    displacements = 'disp_x disp_y disp_z'
    rotations = 'rot_x rot_y rot_z'
    thickness = 0.01
    through_thickness_order = SECOND
  [../]
  [./stress]
    type = ADComputeShellStress2
    block = '100'
    through_thickness_order = SECOND
  [../]
[]

[Postprocessors]
  [./disp_z]
    type = PointValue
    point = '0 0 1'
    variable = disp_z
  [../]
  [./disp_y]
    type = PointValue
    point = '0 1 0'
    variable = disp_y
  [../]
[]

[Outputs]
  exodus = true
[]
