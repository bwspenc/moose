[GlobalParams]
  order = SECOND
  family = RATIONAL_BERNSTEIN
  displacements = 'disp_x disp_y disp_z'
[]

#[Problem]
#  kernel_coverage_check = false
#  material_coverage_check = false
#[]

[Mesh]
  [igafile]
    type = FileMeshGenerator
    file = rpv_simple_uspline_named_sidesets.e
    clear_spline_nodes = true
  []
  [pin]
    type = ExtraNodesetGenerator
    input = igafile
    new_boundary = pin
    nodes = '0'
  []
[]

[Variables]
  [temp]
    initial_condition = 560
  []
[]

[AuxVariables]
  [hoop_stress_clad]
    order = SECOND
    family = MONOMIAL
  []
  [axial_stress_clad]
    order = SECOND
    family = MONOMIAL
  []
  [hoop_stress_base]
    order = SECOND
    family = MONOMIAL
  []
  [axial_stress_base]
    order = SECOND
    family = MONOMIAL
  []
[]

[Modules/TensorMechanics/Master]
  [all]
    strain = SMALL
    temperature = temp
    add_variables = true
    generate_output = 'stress_xx stress_yy stress_zz'
    eigenstrain_names = thermal_eigenstrain
    block = '0 1'
  []
[]

[Kernels]
  [heat]
    type = HeatConduction
    variable = temp
    block = '0 1'
  []
  [heat_ie]
    type = HeatConductionTimeDerivative
    variable = temp
    block = '0 1'
  []
[]

[AuxKernels]
  [axial_stress_clad]
    type = MaterialRealAux
    block = 1
    property = axial_stress_clad
    variable = axial_stress_clad
  []
  [hoop_stress_clad]
    type = MaterialRealAux
    block = 1
    property = hoop_stress_clad
    variable = hoop_stress_clad
  []
  [axial_stress_base]
    type = MaterialRealAux
    block = 0
    property = axial_stress_base
    variable = axial_stress_base
  []
  [hoop_stress_base]
    type = MaterialRealAux
    block = 0
    property = hoop_stress_base
    variable = hoop_stress_base
  []
[]

[Functions]
  [timeFunc]
    type = PiecewiseLinear
    xy_data = '0 1
               1 59
               60 60'
  []
  [coolant_pressure_history]
    type = PiecewiseLinear
    x = '0 500 1200 1220 1700 2000'
    y = '15 3   2  9   2    1'
  []
  [coolant_temperature_history]
    type = PiecewiseLinear
    x = '0 180 2000'
    y = '560 325 425'
  []
[]

[BCs]
  [no_x_all]
    type = DirichletBC
    variable = disp_x
    boundary = 'X0_Plane'
    value = 0.0
  []

  [no_y_all]
    type = DirichletBC
    variable = disp_y
    boundary = 'Y0_Plane'
    value = 0.0
  []
  [pin_z]
    type = DirichletBC
    variable = disp_z
    boundary = 'pin'
    value = 0.0
  []

  [Pressure]
    [coolantPressure]
      boundary = 'Inner_Liner'
      factor = 1e6
      function = coolant_pressure_history
    []
  []
#  [high_T_inner_surface]
#    type = DirichletBC
#    variable = temp
#    boundary = 'Inner_Liner'
#    value = 5000
#  []
  [convective_flux_inner_surface]
    type = ConvectiveFluxFunction
    boundary = 'Inner_Liner'
    variable = temp
    coefficient = 5000.0
    T_infinity = coolant_temperature_history
  []
[]

[Materials]
  [thermal_base]
    type = HeatConductionMaterial
    block = '0'
    thermal_conductivity_temperature_function = 40
    specific_heat_temperature_function =  500
    temp = temp
  []

  [thermal_clad]
    type = HeatConductionMaterial
    block = '1'
    thermal_conductivity_temperature_function = 20
    specific_heat_temperature_function = 500
    temp = temp
  []

  [youngs_modulus_base]
    type = PiecewiseLinearInterpolationMaterial
    x = '200 700'
    y = '29000 29000'
    scale_factor = 6894757.3
    property = youngs_modulus_base
    variable = temp
    block = 0
  []
  [elastic_tensor_base]
    type = ComputeVariableIsotropicElasticityTensor
    args = temp
    youngs_modulus = youngs_modulus_base
    poissons_ratio = 0.3
    block = 0
  []
  [stress_base]
    type = ComputeLinearElasticStress
    block = 0
  []
  [thermal_strain_base]
    type = ComputeMeanThermalExpansionFunctionEigenstrain
    temperature = temp
    thermal_expansion_function = 0.000006
    thermal_expansion_function_reference_temperature = 290
    stress_free_temperature = 525
    eigenstrain_name = thermal_eigenstrain
    block = 0
  []

  [youngs_modulus_clad]
    type = PiecewiseLinearInterpolationMaterial
    x = '200 700'
    y = '22000 22000'
    scale_factor = 6894757.3
    property = youngs_modulus_clad
    variable = temp
    block = 1
  []
  [elastic_tensor_clad]
    type = ComputeVariableIsotropicElasticityTensor
    args = temp
    youngs_modulus = youngs_modulus_clad
    poissons_ratio = 0.3
    block = 1
  []
  [stress_clad]
    type = ComputeLinearElasticStress
    block = 1
  []
  [thermal_strain_clad]
    type = ComputeMeanThermalExpansionFunctionEigenstrain
    temperature = temp
    thermal_expansion_function = 0.000008
    thermal_expansion_function_reference_temperature = 290
    stress_free_temperature = 525
    eigenstrain_name = thermal_eigenstrain
    block = 1
  []

  [density]
    type = Density
    block = '0 1'
    density = 7800
  []

  [axial_stress_clad]
    type = RankTwoCylindricalComponent
    block = 1
    rank_two_tensor = stress
    property_name = axial_stress_clad
    cylindrical_component = AxialStress
    cylindrical_axis_point1 = '0. 0. 0.'
    cylindrical_axis_point2 = '0. 0. 1.'
  []
  [hoop_stress_clad]
    type = RankTwoCylindricalComponent
    block = 1
    rank_two_tensor = stress
    cylindrical_component = HoopStress
    property_name = hoop_stress_clad
    cylindrical_axis_point1 = '0. 0. 0.'
    cylindrical_axis_point2 = '0. 0. 1.'
  []
  [axial_stress_base]
    type = RankTwoCylindricalComponent
    block = 0
    rank_two_tensor = stress
    cylindrical_component = AxialStress
    property_name = axial_stress_base
    cylindrical_axis_point1 = '0. 0. 0.'
    cylindrical_axis_point2 = '0. 0. 1.'
  []
  [hoop_stress_base]
    type = RankTwoCylindricalComponent
    block = 0
    rank_two_tensor = stress
    cylindrical_component = HoopStress
    property_name = hoop_stress_base
    cylindrical_axis_point1 = '0. 0. 0.'
    cylindrical_axis_point2 = '0. 0. 1.'
  []
[]

[VectorPostprocessors]
  [axial_stress_base_0_0]
    type = LineValueSampler
    num_points = 20
    outputs = vpp
    sort_by = id
    start_point = '2.20332 0 3'
    end_point = '2.41415 0 3'
    variable = axial_stress_base
  []
  [axial_stress_clad_0_0]
    type = LineValueSampler
    num_points = 6
    outputs = vpp
    sort_by = id
    start_point = '2.19714 0 3'
    end_point = '2.20112 0 3'
    variable = axial_stress_clad
  []
  [hoop_stress_base_0_0]
    type = LineValueSampler
    num_points = 20
    outputs = vpp
    sort_by = id
    start_point = '2.20332 0 3'
    end_point = '2.41415 0 3'
    variable = hoop_stress_base
  []
  [hoop_stress_clad_0_0]
    type = LineValueSampler
    num_points = 6
    outputs = vpp
    sort_by = id
    start_point = '2.19714 0 3'
    end_point = '2.20112 0 3'
    variable = hoop_stress_clad
  []
  [temperature_base_0_0]
    type = LineValueSampler
    num_points = 20
    outputs = vpp
    sort_by = id
    start_point = '2.20332 0 3'
    end_point = '2.41415 0 3'
    variable = temp
  []
  [coefs_axial_base]
    type = LeastSquaresFitHistory
    order = 4
    vectorpostprocessor = axial_stress_base_0_0
    x_name = id
    x_scale = 4.562
    x_shift = 0.00621538
    y_name = axial_stress_base
  []
  [coefs_axial_clad]
    type = LeastSquaresFitHistory
    order = 1
    vectorpostprocessor = axial_stress_clad_0_0
    x_name = id
    x_scale = 4.562
    x_shift = 4.064e-05
    y_name = axial_stress_clad
  []
  [coefs_hoop_base]
    type = LeastSquaresFitHistory
    order = 4
    vectorpostprocessor = hoop_stress_base_0_0
    x_name = id
    x_scale = 4.562
    x_shift = 0.00621538
    y_name = hoop_stress_base
  []
  [coefs_hoop_clad]
    type = LeastSquaresFitHistory
    order = 1
    vectorpostprocessor = hoop_stress_clad_0_0
    x_name = id
    x_scale = 4.562
    x_shift = 4.064e-05
    y_name = hoop_stress_clad
  []
  [coefs_temp_base]
    type = LeastSquaresFitHistory
    order = 4
    vectorpostprocessor = temperature_base_0_0
    x_name = id
    x_scale = 4.562
    x_shift = 0.00621538
    y_name = temp
  []
[]

[Preconditioning]
  [SMP]
    type = SMP
    full = true
  []
[]

[Executioner]
  automatic_scaling = true
  solve_type = 'PJFNK'
  type = Transient
  petsc_options = '-ksp_snes_ew'

  # Works great on small problems, takes too much RAM on large
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu       superlu_dist'

  # Lousy convergence even on small problems
  # petsc_options_iname = '-ksp_gmres_restart -pc_type -pc_hypre_type -pc_hypre_boomeramg_max_iter'
  # petsc_options_value = ' 201                hypre    boomeramg      4'

  # Slow even on small problems
  # petsc_options_iname = '-ksp_gmres_restart -pc_type -pc_asm_overlap -sub_pc_type -sub_pc_factor_levels'
  # petsc_options_value = ' 201                asm      4               ilu          4'

  l_max_its = 25
  nl_max_its = 50
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-11

  start_time = 0.0
  dt = 1
  end_time = 2000

  [Predictor]
    type = SimplePredictor
    scale = 1.0
    skip_times_old = '1'
  []

  dtmax = 60
  dtmin = 1

  [TimeStepper]
    type = FunctionDT
    function = timeFunc
  []
[]

[Dampers]
  [limitT]
    type = MaxIncrement
    max_increment = 50.0
    variable = temp
  []
[]

[Outputs]
  exodus = true
  csv = true
  [vpp]
    type = CSV
    file_base = 'csv/out'
    execute_on = timestep_end
  []
[]
