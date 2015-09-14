[Mesh]
  file = sliding_elastic_blocks_2d.e
  displacements = 'disp_x disp_y'
[]

#[Problem]
#  type = FrictionalContactProblem
#  master = '2'
#  slave = '3'
#  friction_coefficient = '0.25'
#  slip_factor = '1.0'
#  slip_too_far_factor = '1.0'
#  disp_x = disp_x
#  disp_y = disp_y
#  residual_x = saved_x
#  residual_y = saved_y
#  diag_stiff_x = diag_saved_x
#  diag_stiff_y = diag_saved_y
#  inc_slip_x = inc_slip_x
#  inc_slip_y = inc_slip_y
#  contact_slip_tolerance_factor = 100
#  target_relative_contact_residual = 1.e-4
#  maximum_slip_iterations = 500
#  minimum_slip_iterations = 1
#  slip_updates_per_iteration = 5
#  solution_variables = 'disp_x disp_y'
#  reference_residual_variables = 'saved_x saved_y'
#[]

[Variables]
  [./disp_x]
  [../]
  [./disp_y]
  [../]
[]

[AuxVariables]
  [./penetration]
    order = FIRST
    family = LAGRANGE
  [../]
  [./saved_x]
  [../]
  [./saved_y]
  [../]
  [./diag_saved_x]
  [../]
  [./diag_saved_y]
  [../]
  [./inc_slip_x]
  [../]
  [./inc_slip_y]
  [../]
  [./accum_slip]
  [../]
  [./tang_force_x]
  [../]
  [./tang_force_y]
  [../]
[]

[Functions]
  [./vertical_movement]
    type = ParsedFunction
    value = -t
  [../]
[]

[SolidMechanics]
  [./solid]
    disp_x = disp_x
    disp_y = disp_y
    save_in_disp_y = saved_y
    save_in_disp_x = saved_x
    diag_save_in_disp_y = diag_saved_y
    diag_save_in_disp_x = diag_saved_x
  [../]
[]

[AuxKernels]
  [./inc_slip_x]
    type = PenetrationAux
    variable = inc_slip_x
    execute_on = timestep_end
    quantity = incremental_slip_x
    boundary = 3
    paired_boundary = 2
  [../]
  [./inc_slip_y]
    type = PenetrationAux
    variable = inc_slip_y
    execute_on = timestep_end
    quantity = incremental_slip_y
    boundary = 3
    paired_boundary = 2
  [../]
  [./accum_slip]
    type = PenetrationAux
    variable = accum_slip
    execute_on = timestep_end
    quantity = accumulated_slip
    boundary = 3
    paired_boundary = 2
  [../]
  [./tangential_force_x]
    type = PenetrationAux
    variable = tang_force_x
    execute_on = timestep_end
    quantity = tangential_force_x
    boundary = 3
    paired_boundary = 2
  [../]
  [./tangential_force_y]
    type = PenetrationAux
    variable = tang_force_y
    execute_on = timestep_end
    quantity = tangential_force_y
    boundary = 3
    paired_boundary = 2
  [../]
  [./penetration]
    type = PenetrationAux
    variable = penetration
    boundary = 3
    paired_boundary = 2
  [../]
[]

[Postprocessors]
  [./bot_react_x]
    type = NodalSum
    variable = saved_x
    boundary = 1
  [../]
  [./bot_react_y]
    type = NodalSum
    variable = saved_y
    boundary = 1
  [../]
  [./top_react_x]
    type = NodalSum
    variable = saved_x
    boundary = 4
  [../]
  [./top_react_y]
    type = NodalSum
    variable = saved_y
    boundary = 4
  [../]
  [./ref_resid_x]
    type = NodalL2Norm
    execute_on = timestep_end
    variable = saved_x
  [../]
  [./ref_resid_y]
    type = NodalL2Norm
    execute_on = timestep_end
    variable = saved_y
  [../]
[]

[BCs]
  [./left_x]
    type = DirichletBC
    variable = disp_x
    boundary = 1
    value = 0.0
  [../]
  [./left_y]
    type = DirichletBC
    variable = disp_y
    boundary = 1
    value = 0.0
  [../]
  [./right_x]
    type = PresetBC
    variable = disp_x
    boundary = 4
    value = -0.005
  [../]
  [./right_y]
    type = FunctionPresetBC
    variable = disp_y
    boundary = 4
    function = vertical_movement
  [../]
[]

[Materials]
  [./left]
    type = LinearIsotropicMaterial
    block = 1
    disp_y = disp_y
    disp_x = disp_x
    poissons_ratio = 0.3
    youngs_modulus = 1e7
  [../]
  [./right]
    type = LinearIsotropicMaterial
    block = 2
    disp_y = disp_y
    disp_x = disp_x
    poissons_ratio = 0.3
    youngs_modulus = 1e6
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

  #Preconditioned JFNK (default)
  solve_type = 'PJFNK'



  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu     superlu_dist'


  line_search = 'none'
#  line_search = 'cp'

  nl_abs_tol = 1e-7
  l_max_its = 100
  nl_max_its = 1000
  dt = 0.01
  end_time = 0.05
  num_steps = 1000
  nl_rel_tol = 1e-6
  dtmin = 0.01
  l_tol = 1e-3

#  [./Predictor]
#    type = SimplePredictor
#    scale = 1.0
#  [../]
[]

[Outputs]
  file_base = slid_elas_blk_2d_new_kin_out
  output_initial = true
  print_linear_residuals = true
  print_perf_log = true
  [./exodus]
    type = Exodus
    elemental_as_nodal = true
  [../]
  [./console]
    type = Console
    max_rows = 5
  [../]
[]

[Contact]
  [./leftright]
    slave = 3
    disp_y = disp_y
    disp_x = disp_x
    master = 2
    model = coulomb
    system = constraint
    friction_coefficient = '0.25'
    penalty = 2e+5
  [../]
[]

[Problem]
#  type = ReferenceResidualProblem
  type = FrictionalContactDamperProblem
  master = '2'
  slave = '3'
  disp_x = disp_x
  disp_y = disp_y
  solution_variables = 'disp_x disp_y'
  reference_residual_variables = 'saved_x saved_y'
[]
