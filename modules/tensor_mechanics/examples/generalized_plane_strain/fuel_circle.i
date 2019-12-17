[GlobalParams]
  displacements = 'disp_x disp_y'
  scalar_out_of_plane_strain = nl_strain_zz
[]

[Problem]
  type = ReferenceResidualProblem
  extra_tag_vectors = 'ref'
  reference_vector = 'ref'
[]

[Mesh]
  [./fmg]
    type = FileMeshGenerator
    file = fuel_circle.e
  [../]
[]

[Variables]
  [./disp_x]
  [../]
  [./disp_y]
  [../]
  [./nl_strain_zz]
    family = SCALAR
    order = FIRST
  [../]
[]

[AuxVariables]
  [./temp]
  [../]
  [./hoop_stress]
#    family = MONOMIAL
#    order = FIRST
  [../]
  [./axial_stress]
#    family = MONOMIAL
#    order = FIRST
  [../]
[]

[BCs]
  [./center_x]
    type = PresetBC
    variable = disp_x
    boundary = 1
    value = 0.0
  [../]
  [./center_y]
    type = PresetBC
    variable = disp_y
    boundary = 1
    value = 0.0
  [../]
  [./right_y]
    type = PresetBC
    variable = disp_y
    boundary = 2
    value = 0.0
  [../]
[]

[Modules/TensorMechanics/Master]
  [gps]
    planar_formulation = GENERALIZED_PLANE_STRAIN
    strain = SMALL
    generate_output = 'stress_xx stress_yy stress_zz strain_xx strain_yy strain_zz vonmises_stress'
    extra_vector_tags = 'ref'
    eigenstrain_names = 'thermal_eigenstrain'
    temperature = temp
  []
[]

[AuxKernels]
  [./temp]
    type = FunctionAux
    variable = temp
    function = 'r:=sqrt(x^2+y^2);r^2*(800-1200)/(0.004^2)+1200'
  [../]
  [./hoop_stress]
    type = RankTwoScalarAux
    variable = hoop_stress
    scalar_type = HoopStress
    rank_two_tensor = stress
    point1 = '0 0 0'
    point2 = '0 0 1'
    execute_on = timestep_end
  [../]
  [./axial_stress]
    type = RankTwoScalarAux
    variable = axial_stress
    scalar_type = AxialStress
    rank_two_tensor = stress
    point1 = '0 0 0'
    point2 = '0 0 1'
    execute_on = timestep_end
  [../]
[]

[Materials]
  [./elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 1e9
    poissons_ratio = 0.3
  [../]
  [./stress]
    type = ComputeLinearElasticStress
  [../]
  [./thermal_strain]
    type = ComputeThermalExpansionEigenstrain
    temperature = temp
    thermal_expansion_coeff = 1e-5
    stress_free_temperature = 293
    eigenstrain_name = thermal_eigenstrain
  [../]
[]

[VectorPostprocessors]
  [disp]
    type = LineValueSampler
    variable = 'disp_x'
    start_point = '0 0 0'
    end_point = '0.004 0 0'
    num_points = 20
    sort_by = id
  []
  [hoop_stress]
    type = LineValueSampler
    variable = 'hoop_stress'
    start_point = '0 0 0'
    end_point = '0.004 0 0'
    num_points = 20
    sort_by = id
    execute_on = timestep_end
  []
  [axial_stress]
    type = LineValueSampler
    variable = 'hoop_stress'
    start_point = '0 0 0'
    end_point = '0.004 0 0'
    num_points = 20
    sort_by = id
    execute_on = timestep_end
  []
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
  nl_abs_tol = 1e-9

  petsc_options = '-snes_ksp_ew'
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu superlu_dist'

  end_time = 1.0
#  automatic_scaling = true #TODO: Doesn't work with GPS!
[]

[Outputs]
  exodus = true
  csv = true
[]
