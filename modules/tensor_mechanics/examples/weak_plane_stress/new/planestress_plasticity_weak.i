[GlobalParams]
  displacements = 'disp_x disp_y'
  out_of_plane_strain = nl_strain_zz
  # volumetric_locking_correction = true
[]

[Problem]
  type = ReferenceResidualProblem
  extra_tag_vectors = 'ref'
  reference_vector = 'ref'
[]

[Mesh]
  [./fmg]
    type = FileMeshGenerator
    file = platehole_04_q4.e
  [../]
[]

[Variables]
  [./disp_x]
  [../]
  [./disp_y]
  [../]
  [./nl_strain_zz]
#    family = MONOMIAL
#    order = CONSTANT
  [../]
[]

[AuxVariables]
  [./plastic_strain_eff]
    order = CONSTANT
    family = MONOMIAL
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
  [./right_x]
    type = FunctionDirichletBC
    variable = disp_x
    boundary = 2
    function = '0.01*t'
  [../]
[]

[Modules/TensorMechanics/Master]
  [plane_stress]
    planar_formulation = WEAK_PLANE_STRESS
    strain = FINITE
    generate_output = 'stress_xx stress_yy stress_zz strain_xx strain_yy strain_zz vonmises_stress'
    decomposition_method = EigenSolution
    extra_vector_tags = 'ref'
  []
[]

[AuxKernels]
  [./plastic_strain_eff]
    type = MaterialRealAux
    property = effective_plastic_strain
    variable = plastic_strain_eff
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
    hardening_constant = 1e4
  [../]
[]

[Postprocessors]
  [./react_z]
    type = MaterialTensorIntegral
    rank_two_tensor = stress
    index_i = 2
    index_j = 2
  [../]
  [./min_strain_zz]
    type = ElementExtremeValue
    variable = strain_zz
    value_type = min
  [../]
  [./max_strain_zz]
    type = ElementExtremeValue
    variable = strain_zz
    value_type = max
  [../]
  [./min_stress_zz]
    type = ElementExtremeValue
    variable = stress_zz
    value_type = min
  [../]
  [./max_stress_zz]
    type = ElementExtremeValue
    variable = stress_zz
    value_type = max
  [../]
  [./max_vonmises_stress]
    type = ElementExtremeValue
    variable = vonmises_stress
    value_type = max
  [../]
  [./ave_vonmises_stress]
    type = ElementAverageValue
    variable = stress_xx
  [../]
  [./stress_zz_diff]
    type = RankTwoComponentL2Difference
    rank_two_tensor = stress
    function = 0
    index_i = 2
    index_j = 2
  [../]
  [./stress_zz_rel_err]
    type = RatioPostprocessor
    dividend = stress_zz_diff
    divisor = ave_vonmises_stress
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
  csv = true
[]
