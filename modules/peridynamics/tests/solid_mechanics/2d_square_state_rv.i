# Test for ordinary state-based peridynamic formulation
# for regular grid from generated mesh with varying bond constants
# partial Jacobian
# Jacobian from bond-based formulation is used for preconditioning

# Square plate with Dirichlet boundary conditions applied
# at the left, top and bottom edges

[GlobalParams]
  displacements = 'disp_x disp_y'
  bond_status = bond_status
  bond_contact = bond_contact
  bond_contact_strain = bond_contact_strain
[]

[Mesh]
  type = GeneratedMeshPD
  horizon_number = 3
  pddim = 2
  shape = 1               # 1. Rectangular, 0. Disk.
  nx = 10
  partitioner = linear
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
  [./bond_contact]
    order = CONSTANT
    family = MONOMIAL
    initial_condition = 0.0
  [../]
  [./bond_critical_strain]
    order = CONSTANT
    family = MONOMIAL
    initial_condition = 1.0
  [../]
  [./bond_contact_strain]
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
  [./Peridynamics]
    state_based_formulation = true
  [../]
[]

[Materials]
  [./linelast]
    type = VElasticSPDMaterial
    youngs_modulus = 2e5
    poissons_ratio = 0.0
  [../]
[]

[Outputs]
  file_base = 2d_square_state_rv
  exodus = true
  print_linear_residuals = true
[]
