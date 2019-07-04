[GlobalParams]
  order = FIRST
  family = LAGRANGE
  displacements = 'disp_x'
[]

[Problem]
  coord_type = RZ
[]

[MeshGenerators]
  [gmg]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 20
    xmin = 2.2
    xmax = 2.4
  []
[]

[Mesh]
  type = MeshGeneratorMesh
[]

[Variables]
  [./disp_x]
  [../]
  [./temp]
    initial_condition = 400
  [../]
[]

[AuxVariables]
  [./hoop_stress]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[Functions]
  [./temp_inner_func]
    type = PiecewiseLinear
    xy_data = '0 400
               1 350'
  [../]
  [./temp_outer_func]
    type = PiecewiseLinear
    xy_data = '0 400
               1 400'
  [../]
  [./press_func]
    type = PiecewiseLinear
    xy_data = '0 15
               1 15'
  [../]
[]

[Kernels]
  [./TensorMechanics]
    use_displaced_mesh = true
    temperature = temp
    #generate_output = 'stress_xx stress_yy stress_zz vonmises_stress hydrostatic_stress'
  [../]
  [./heat]
    type = HeatConduction
    variable = temp
  [../]
[]

#[Modules/TensorMechanics/Master]
#  [./all]
#    incremental = true
#    strain = FINITE
#    eigenstrain_names = thermal_expansion
#    generate_output = 'stress_xx stress_yy stress_zz vonmises_stress hydrostatic_stress'
#  [../]
#[]

[AuxKernels]
  [./hoop_stress]
    type = RankTwoScalarAux
    rank_two_tensor = stress
    variable = hoop_stress
    scalar_type = HoopStress
    execute_on = timestep_end
  [../]
[]

[BCs]
  [./Pressure]
    [./internal_pressure]
      boundary = 'left'
      factor = 1.e6
      function = press_func
    [../]
  [../]

  [./t_in]
    type = FunctionDirichletBC
    variable = temp
    boundary = 'left'
    function = temp_inner_func
  [../]

  [./t_out]
    type = FunctionDirichletBC
    variable = temp
    boundary = 'right'
    function = temp_outer_func
  [../]
[]

[Materials]
  [./thermal1]
    type = HeatConductionMaterial
    thermal_conductivity = 25.0
    specific_heat = 490.0
    temp = temp
  [../]

  [./elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 193.05e9
    poissons_ratio = 0.3
  [../]

  [./strain]
    type = ComputeAxisymmetric1DFiniteStrain
    eigenstrain_names = thermal_expansion
    out_of_plane_strain = 0
  [../]

  [./stress]
    type = ComputeFiniteStrainElasticStress
  [../]

  [./thermal_expansion]
    type = ComputeThermalExpansionEigenstrain
    thermal_expansion_coeff = 13e-6
    stress_free_temperature = 295.00
    temperature = temp
    eigenstrain_name = thermal_expansion
  [../]

  [./density]
    type = Density
    density = 8000.0
  [../]
[]

[Executioner]
  type = Transient
  solve_type = 'PJFNK'

  l_max_its = 25
  nl_max_its = 20
  nl_rel_tol = 1e-9
  l_tol = 1e-2

  start_time = 0.0
  dt = 1
  end_time = 1
  dtmin = 1
[]

[Outputs]
  file_base = 1d_out
  exodus = true
[]
