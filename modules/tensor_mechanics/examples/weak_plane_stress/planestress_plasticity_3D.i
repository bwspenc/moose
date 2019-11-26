[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
  # volumetric_locking_correction = true
[]

[Problem]
  type = ReferenceResidualProblem
  extra_tag_vectors = 'ref'
  reference_vector = 'ref'
[]

[Mesh]
  file = 3D_hex20_3470.e
[]

[Variables]
  [./disp_x]
    order = SECOND
  [../]
  [./disp_y]
    order = SECOND
  [../]
  [./disp_z]
    order = SECOND
  [../]
[]

[AuxVariables]
  [./strain_xx]
    family = MONOMIAL
    order = SECOND
  [../]
  [./strain_yy]
    family = MONOMIAL
    order = SECOND
  [../]
  [./strain_zz]
    family = MONOMIAL
    order = SECOND
  [../]
  [./total_strain_zz]
    family = MONOMIAL
    order = SECOND
  [../]
  [./strain_xy]
    family = MONOMIAL
    order = SECOND
  [../]

  [./plastic_strain_eff]
    family = MONOMIAL
    order = CONSTANT
  [../]

  [./stress_xx]
    family = MONOMIAL
    order = SECOND
  [../]
  [./stress_yy]
    family = MONOMIAL
    order = SECOND
  [../]
  [./stress_xy]
    family = MONOMIAL
    order = SECOND
  [../]
  [./stress_zz]
    family = MONOMIAL
    order = SECOND
  [../]
  [./vonmises]
    family = MONOMIAL
    order = SECOND
  [../]
[]

[BCs]
  [./left_x]
    type = PresetBC
    variable = disp_x
    boundary = 4
    value = 0.0
  [../]
  [./bottom_y]
    type = PresetBC
    variable = disp_y
    boundary = 1
    value = 0.0
  [../]
  [./back_z]
    type = PresetBC
    variable = disp_z
    boundary = 5
    value = 0.0
  [../]
  # [./fix_z]
  #   type = PresetBC
  #   variable = disp_z
  #   boundary = 6
  #   value = 0.0
  # [../]
  [./right_x]
    type = FunctionDirichletBC
    variable = disp_x
    boundary = 2
    function = '0.01*t'
  [../]
  # [./top_y]
  #   type = FunctionDirichletBC
  #   variable = disp_y
  #   boundary = 3
  #   function = '0.01*t'
  # [../]
[]

[Modules/TensorMechanics/Master]
  [./all]
    strain = FINITE
    decomposition_method = EigenSolution
    extra_vector_tags = 'ref'
  [../]
[]

[AuxKernels]
  [./strain_xx]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    variable = strain_xx
    index_i = 0
    index_j = 0
  [../]
  [./strain_yy]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    variable = strain_yy
    index_i = 1
    index_j = 1
  [../]
  [./strain_zz]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    variable = strain_zz
    index_i = 2
    index_j = 2
  [../]
  [./total_strain_zz]
    type = RankTwoAux
    rank_two_tensor = total_strain
    variable = total_strain_zz
    index_i = 2
    index_j = 2
  [../]
  [./strain_xy]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    variable = strain_xy
    index_i = 0
    index_j = 1
  [../]

  [./plastic_strain_eff]
    type = MaterialRealAux
    property = effective_plastic_strain
    variable = plastic_strain_eff
  [../]

  [./stress_xx]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_xx
    index_i = 0
    index_j = 0
  [../]
  [./stress_yy]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_yy
    index_i = 1
    index_j = 1
  [../]
  [./stress_xy]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_xy
    index_i = 0
    index_j = 1
  [../]
  [./stress_zz]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_zz
    index_i = 2
    index_j = 2
  [../]
  [./vonmises]
    type = RankTwoScalarAux
    rank_two_tensor = stress
    variable = vonmises
    scalar_type = VonMisesStress
    execute_on = timestep_end
  [../]
[]

[Materials]
  [./elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 1e6
    poissons_ratio = 0.3
  [../]

  [./stress]
    type = ComputeMultipleInelasticStress
    inelastic_models = 'isoplasticity'
  [../]
  [./isoplasticity]
    type = IsotropicPlasticityStressUpdate
    yield_stress = 2e3
    hardening_constant = 1e3
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
  nl_max_its = 100

  petsc_options = '-snes_ksp_ew'
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu superlu_dist'
  line_search = none

  start_time = 0
  dt = 0.02
  dtmin = 0.01
  end_time = 1.0
  automatic_scaling = true
  [Predictor]
    type = SimplePredictor
    scale = 1.0
  []
[]

[Outputs]
  exodus = true
[]
