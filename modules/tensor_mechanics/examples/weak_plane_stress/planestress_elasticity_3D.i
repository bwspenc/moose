
[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Mesh]
  file = 3D_3470.e
[]

[Variables]
  [./disp_x]
  [../]
  [./disp_y]
  [../]
  [./disp_z]
  [../]
[]

[AuxVariables]
  [./strain_xx]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./strain_yy]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./strain_zz]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./strain_xy]
    family = MONOMIAL
    order = CONSTANT
  [../]

  [./stress_xx]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./stress_yy]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./stress_xy]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./stress_zz]
    family = MONOMIAL
    order = CONSTANT
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
  [./right_x]
    type = FunctionPresetBC
    variable = disp_x
    boundary = 2
    function = '7.280419e-4*(x+0.5)/((x+0.5)^2+(y+0.5)^2)*(13.735473*((x+0.5)^2+(y+0.5)^2)+6.153846+(2-0.448026/((x+0.5)^2+(y+0.5)^2))*(4*(x+0.5)*(x+0.5)/((x+0.5)^2+(y+0.5)^2)-3))'
  [../]
  [./right_y]
    type = FunctionPresetBC
    variable = disp_y
    boundary = 2
    function = '7.280419e-4*(y+0.5)/((x+0.5)^2+(y+0.5)^2)*(-4.120642*((x+0.5)^2+(y+0.5)^2)-2.153846+(2-0.448026/((x+0.5)^2+(y+0.5)^2))*(3-4*(y+0.5)*(y+0.5)/((x+0.5)^2+(y+0.5)^2)))'
  [../]
  [./top_x]
    type = FunctionPresetBC
    variable = disp_x
    boundary = 3
    function = '7.280419e-4*(x+0.5)/((x+0.5)^2+(y+0.5)^2)*(13.735473*((x+0.5)^2+(y+0.5)^2)+6.153846+(2-0.448026/((x+0.5)^2+(y+0.5)^2))*(4*(x+0.5)*(x+0.5)/((x+0.5)^2+(y+0.5)^2)-3))'
  [../]
  [./top_y]
    type = FunctionPresetBC
    variable = disp_y
    boundary = 3
    function = '7.280419e-4*(y+0.5)/((x+0.5)^2+(y+0.5)^2)*(-4.120642*((x+0.5)^2+(y+0.5)^2)-2.153846+(2-0.448026/((x+0.5)^2+(y+0.5)^2))*(3-4*(y+0.5)*(y+0.5)/((x+0.5)^2+(y+0.5)^2)))'
  [../]
[]

[Modules/TensorMechanics/Master]
  [./all]
  [../]
[]

[AuxKernels]
  [./strain_xx]
    type = RankTwoAux
    rank_two_tensor = total_strain
    variable = strain_xx
    index_i = 0
    index_j = 0
  [../]
  [./strain_yy]
    type = RankTwoAux
    rank_two_tensor = total_strain
    variable = strain_yy
    index_i = 1
    index_j = 1
  [../]
  [./strain_zz]
    type = RankTwoAux
    rank_two_tensor = total_strain
    variable = strain_zz
    index_i = 2
    index_j = 2
  [../]
  [./strain_xy]
    type = RankTwoAux
    rank_two_tensor = total_strain
    variable = strain_xy
    index_i = 0
    index_j = 1
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
[]

[Materials]
  [./elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 1e6
    poissons_ratio = 0.3
  [../]
  [./stress]
    type = ComputeLinearElasticStress
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

  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu superlu_dist'

  end_time = 1.0
[]

[Outputs]
  exodus = true
[]
