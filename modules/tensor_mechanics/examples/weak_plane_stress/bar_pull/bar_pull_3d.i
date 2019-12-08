[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Problem]
  type = ReferenceResidualProblem
  extra_tag_vectors = 'ref'
  reference_vector = 'ref'
[]

[Mesh]
  [./gmg]
    type = GeneratedMeshGenerator
    dim = 3
    xmin = 0
    ymin = 0
    zmin = 0
    xmax = 5
    ymax = 1
    zmax = 1
    nx = 5
    ny = 1
    nz = 1
  [../]
[]

[Variables]
  [./disp_x]
  [../]
  [./disp_y]
  [../]
  [./disp_z]
  [../]
[]

[BCs]
  [./right_x]
    type = FunctionPresetBC
    variable = disp_x
    boundary = right
    function = t
  [../]
  [./left_x]
    type = PresetBC
    variable = disp_x
    boundary = left
    value = 0.0
  [../]
  [./left_y]
    type = PresetBC
    variable = disp_y
    boundary = left
    value = 0.0
  [../]
  [./left_z]
    type = PresetBC
    variable = disp_z
    boundary = left
    value = 0.0
  [../]
[]

[Modules/TensorMechanics/Master]
  [3d]
    strain = FINITE
    generate_output = 'stress_xx stress_yy stress_zz strain_xx strain_yy strain_zz vonmises_stress'
    decomposition_method = EigenSolution
  []
[]

[Materials]
  [./elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 1e6
    poissons_ratio = 0.3
  [../]
  [./stress]
    type = ComputeFiniteStrainElasticStress
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
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-11

  petsc_options = '-snes_ksp_ew'
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu superlu_dist'

  start_time = 0.0
  dt = 0.1
  end_time = 1.0
  automatic_scaling = true
[]

[Outputs]
  exodus = true
  csv = true
[]
