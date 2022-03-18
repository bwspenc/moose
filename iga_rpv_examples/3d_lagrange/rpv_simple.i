[GlobalParams]
  order = SECOND
  family = LAGRANGE
  displacements = 'disp_x disp_y disp_z'
[]

[Mesh]
  file = rpv_simple.e
  use_displaced_mesh = false
[]

[Problem]
  type = ReferenceResidualProblem
  extra_tag_vectors = 'ref'
  reference_vector = 'ref'
  group_variables = 'disp_x disp_y disp_z'
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
    extra_vector_tags = 'ref'
  []
[]

[Kernels]
  [heat]
    type = HeatConduction
    variable = temp
    extra_vector_tags = 'ref'
  []
  [heat_ie]
    type = HeatConductionTimeDerivative
    variable = temp
    extra_vector_tags = 'ref'
  []
[]

[AuxKernels]
  [axial_stress_clad]
    type = MaterialRealAux
    block = 2
    property = axial_stress_clad
    variable = axial_stress_clad
  []
  [hoop_stress_clad]
    type = MaterialRealAux
    block = 2
    property = hoop_stress_clad
    variable = hoop_stress_clad
  []
  [axial_stress_base]
    type = MaterialRealAux
    block = 1
    property = axial_stress_base
    variable = axial_stress_base
  []
  [hoop_stress_base]
    type = MaterialRealAux
    block = 1
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
    boundary = '1'
    value = 0.0
  []

  [no_y_all]
    type = DirichletBC
    variable = disp_y
    boundary = '2'
    value = 0.0
  []
  [no_z]
    type = DirichletBC
    variable = disp_z
    boundary = '101'
    value = 0.0
  []

  [Pressure]
    [coolantPressure]
      boundary = 3
      factor = 1e6
      function = coolant_pressure_history
    []
  []
  [convective_flux_inner_surface]
    type = ConvectiveFluxFunction
    boundary = 3
    variable = temp
    coefficient = 5000.0
    T_infinity = coolant_temperature_history
  []
[]

[Materials]
  [thermal_base]
    type = HeatConductionMaterial
    block = '1'
    thermal_conductivity_temperature_function = 40
    specific_heat_temperature_function = 500
    temp = temp
  []

  [thermal_clad]
    type = HeatConductionMaterial
    block = '2'
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
    block = 1
  []
  [elastic_tensor_base]
    type = ComputeVariableIsotropicElasticityTensor
    args = temp
    youngs_modulus = youngs_modulus_base
    poissons_ratio = 0.3
    block = 1
  []
  [stress_base]
    type = ComputeLinearElasticStress
    block = 1
  []
  [thermal_strain_base]
    type = ComputeMeanThermalExpansionFunctionEigenstrain
    temperature = temp
    thermal_expansion_function = 0.000006
    thermal_expansion_function_reference_temperature = 290
    stress_free_temperature = 525
    eigenstrain_name = thermal_eigenstrain
    block = 1
  []

  [youngs_modulus_clad]
    type = PiecewiseLinearInterpolationMaterial
    x = '200 700'
    y = '22000 22000'
    scale_factor = 6894757.3
    property = youngs_modulus_clad
    variable = temp
    block = 2
  []
  [elastic_tensor_clad]
    type = ComputeVariableIsotropicElasticityTensor
    args = temp
    youngs_modulus = youngs_modulus_clad
    poissons_ratio = 0.3
    block = 2
  []
  [stress_clad]
    type = ComputeLinearElasticStress
    block = 2
  []
  [thermal_strain_clad]
    type = ComputeMeanThermalExpansionFunctionEigenstrain
    temperature = temp
    thermal_expansion_function = 0.000008
    thermal_expansion_function_reference_temperature = 290
    stress_free_temperature = 525
    eigenstrain_name = thermal_eigenstrain
    block = 2
  []

  [density]
    type = Density
    block = '1 2'
    density = 7800.0
  []

  [axial_stress_clad]
    type = RankTwoCylindricalComponent
    block = 2
    rank_two_tensor = stress
    property_name = axial_stress_clad
    cylindrical_component = AxialStress
    cylindrical_axis_point1 = '0. 0. 0.'
    cylindrical_axis_point2 = '0. 0. 1.'
  []
  [hoop_stress_clad]
    type = RankTwoCylindricalComponent
    block = 2
    rank_two_tensor = stress
    cylindrical_component = HoopStress
    property_name = hoop_stress_clad
    cylindrical_axis_point1 = '0. 0. 0.'
    cylindrical_axis_point2 = '0. 0. 1.'
  []
  [axial_stress_base]
    type = RankTwoCylindricalComponent
    block = 1
    rank_two_tensor = stress
    cylindrical_component = AxialStress
    property_name = axial_stress_base
    cylindrical_axis_point1 = '0. 0. 0.'
    cylindrical_axis_point2 = '0. 0. 1.'
  []
  [hoop_stress_base]
    type = RankTwoCylindricalComponent
    block = 1
    rank_two_tensor = stress
    cylindrical_component = HoopStress
    property_name = hoop_stress_base
    cylindrical_axis_point1 = '0. 0. 0.'
    cylindrical_axis_point2 = '0. 0. 1.'
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
  # petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  # petsc_options_value = 'lu       superlu_dist'
  petsc_options_iname = '-ksp_gmres_restart -pc_type -pc_hypre_type -pc_hypre_boomeramg_max_iter'
  petsc_options_value = ' 201                hypre    boomeramg      4'

  l_max_its = 25
  nl_max_its = 50
  nl_rel_tol = 1e-8
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

#[LineSamplers]
#  total_thickness = 0.219202
#  clad_thickness = 0.004064
#  inner_radius = 2.1971
#  base_block = 1
#  clad_block = 2
#  coord_system = 3D_CARTESIAN
#  axis_of_rotation = z
#  line_sampler_type = LineValueSampler
#  temperature = 'temp'
#  axial_stress_base = axial_stress_base
#  hoop_stress_base = hoop_stress_base
#  axial_stress_clad = axial_stress_clad
#  hoop_stress_clad = hoop_stress_clad
#  base_order = 4
#  axial_start = -0.2
#  axial_end = 4.2
#  axial_num_points = 16
#  azimuthal_start = 0
#  azimuthal_end = 90
#  azimuthal_num_points = 12
#  azimuth0 = '1 0 0'
#  base_thickness_offset_fraction = 0.01
#  clad_thickness_offset_fraction = 0.01
#[]

[VectorPostprocessors]
  [axial_stress_base_0_0]
    type = LineValueSampler
    num_points = 20
    outputs = vpp
    sort_by = id
    start_point = '2.20332 0 0'
    end_point = '2.41415 0 0'
    variable = axial_stress_base
  []
  [axial_stress_clad_0_0]
    type = LineValueSampler
    num_points = 6
    outputs = vpp
    sort_by = id
    start_point = '2.19714 0 0'
    end_point = '2.20112 0 0'
    variable = axial_stress_clad
  []
  [hoop_stress_base_0_0]
    type = LineValueSampler
    num_points = 20
    outputs = vpp
    sort_by = id
    start_point = '2.20332 0 0'
    end_point = '2.41415 0 0'
    variable = hoop_stress_base
  []
  [hoop_stress_clad_0_0]
    type = LineValueSampler
    num_points = 6
    outputs = vpp
    sort_by = id
    start_point = '2.19714 0 0'
    end_point = '2.20112 0 0'
    variable = hoop_stress_clad
  []
  [temperature_base_0_0]
    type = LineValueSampler
    num_points = 20
    outputs = vpp
    sort_by = id
    start_point = '2.20332 0 0'
    end_point = '2.41415 0 0'
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

[Dampers]
  [limitT]
    type = MaxIncrement
    max_increment = 50.0
    variable = temp
  []
[]

[Outputs]
  exodus = true
  perf_graph = true
  csv = true
  [vpp]
    type = CSV
    file_base = 'csv/out'
    execute_on = timestep_end
  []
[]
