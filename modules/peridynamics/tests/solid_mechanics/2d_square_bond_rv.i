# Test for bond-based peridynamic formulation
# for regular grid from generated mesh with varying bond constants

# Square plate with Dirichlet boundary conditions applied
# at the left, top and bottom edges

[GlobalParams]
  displacements = 'disp_x disp_y'
  bond_status = bond_status
[]

[Mesh]
  type = GeneratedMeshPD
  horizon_number = 3
  pddim = 2
  shape = 1               # 1. Rectangular, 0. Disk.
  nx = 10
[]

[Variables]
  [./disp_x]
  [../]
  [./disp_y]
  [../]
[]

[AuxVariables]
  [./bond_status]
    order = CONSTANT
    family = MONOMIAL
    initial_condition = 1.0
  [../]
[]

[BCs]
  [./left_dx]
    type = DirichletBC
    variable = disp_x
    boundary = 0
    value = 0.0
  [../]
  [./top_dy]
    type = DirichletBC
    variable = disp_y
    boundary = 3
    value = 0.0
  [../]
  [./bottom_dy]
    type = FunctionPresetBC
    variable = disp_y
    boundary = 2
    function = '-0.001*t'
  [../]
[]

[Preconditioning]
  [./SMP]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient
  solve_type = PJFNK
  start_time = 0
  end_time = 1
[]

[Kernels]
  [./solid_x]
    type = StressDivergenceBPD
    variable = disp_x
    component = 0
  [../]
  [./solid_y]
    type = StressDivergenceBPD
    variable = disp_y
    component = 1
  [../]
[]

[Materials]
  [./linelast]
    type = VElasticBPDMaterial
    youngs_modulus = 2e5
    poissons_ratio = 0.33
  [../]
[]

[Outputs]
  file_base = 2d_square_bond_rv
  exodus = true
  print_linear_residuals = true
[]
