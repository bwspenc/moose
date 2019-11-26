[GlobalParams]
  displacements = 'disp_x disp_y'
  out_of_plane_strain = strain_zz
[]

[Problem]
  type = ReferenceResidualProblem
  extra_tag_vectors = 'ref'
  reference_vector = 'ref'
[]

[Mesh]
  # file = 2D_91.e
  #file = 2D_936.e
  #file = 2D_4901.e
  file = 2D_9353.e
[]

[Variables]
  [./disp_x]
  [../]
  [./disp_y]
  [../]

  [./strain_zz]
    family = MONOMIAL
    order = CONSTANT
  [../]
[]

[AuxVariables]
  [./strain_xx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strain_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strain_xy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./nl_strain_zz]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./stress_xx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_xy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_zz]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[Postprocessors]
  [./react_z]
    type = MaterialTensorIntegral
    rank_two_tensor = stress
    index_i = 2
    index_j = 2
  [../]
#  [./min_strain_zz]
#    type = NodalExtremeValue
#    variable = strain_zz
#    value_type = min
#  [../]
#  [./max_strain_zz]
#    type = NodalExtremeValue
#    variable = strain_zz
#    value_type = max
#  [../]
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
[]

[BCs]
  [./left_x]
    type = PresetBC
    variable = disp_x
    boundary = 1
    value = 0.0
  [../]
  [./bottom_y]
    type = PresetBC
    variable = disp_y
    boundary = 3
    value = 0.0
  [../]
  [./right_x]
    type = FunctionDirichletBC
    variable = disp_x
    boundary = 2
    #function = '7.280419e-4*(x+0.5)/((x+0.5)^2+(y+0.5)^2)*(12.499281*((x+0.5)^2+(y+0.5)^2)+5.6+(2-0.448026/((x+0.5)^2+(y+0.5)^2))*(4*(x+0.5)*(x+0.5)/((x+0.5)^2+(y+0.5)^2)-3))'
    function = '7.280419e-4*(x+0.5)/((x+0.5)^2+(y+0.5)^2)*(13.735473*((x+0.5)^2+(y+0.5)^2)+6.153846+(2-0.448026/((x+0.5)^2+(y+0.5)^2))*(4*(x+0.5)*(x+0.5)/((x+0.5)^2+(y+0.5)^2)-3))'
  [../]
  [./right_y]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = 2
    #function = '7.280419e-4*(y+0.5)/((x+0.5)^2+(y+0.5)^2)*(-5.356835*((x+0.5)^2+(y+0.5)^2)-1.6+(2-0.448026/((x+0.5)^2+(y+0.5)^2))*(3-4*(y+0.5)*(y+0.5)/((x+0.5)^2+(y+0.5)^2)))'
    function = '7.280419e-4*(y+0.5)/((x+0.5)^2+(y+0.5)^2)*(-4.120642*((x+0.5)^2+(y+0.5)^2)-2.153846+(2-0.448026/((x+0.5)^2+(y+0.5)^2))*(3-4*(y+0.5)*(y+0.5)/((x+0.5)^2+(y+0.5)^2)))'
  [../]
  [./top_x]
    type = FunctionDirichletBC
    variable = disp_x
    boundary = 4
    function = '7.280419e-4*(x+0.5)/((x+0.5)^2+(y+0.5)^2)*(13.735473*((x+0.5)^2+(y+0.5)^2)+6.153846+(2-0.448026/((x+0.5)^2+(y+0.5)^2))*(4*(x+0.5)*(x+0.5)/((x+0.5)^2+(y+0.5)^2)-3))'
  [../]
  [./top_y]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = 4
    function = '7.280419e-4*(y+0.5)/((x+0.5)^2+(y+0.5)^2)*(-4.120642*((x+0.5)^2+(y+0.5)^2)-2.153846+(2-0.448026/((x+0.5)^2+(y+0.5)^2))*(3-4*(y+0.5)*(y+0.5)/((x+0.5)^2+(y+0.5)^2)))'
  [../]
[]

[Kernels]
  [./disp_x]
    type = StressDivergenceTensors
    variable = disp_x
    component = 0
    extra_vector_tags = 'ref'
  [../]
  [./disp_y]
    type = StressDivergenceTensors
    variable = disp_y
    component = 1
    extra_vector_tags = 'ref'
  [../]

  [./solid_z]
    type = WeakPlaneStress
    variable = strain_zz
    extra_vector_tags = 'ref'
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
  [./strain_xy]
    type = RankTwoAux
    rank_two_tensor = total_strain
    variable = strain_xy
    index_i = 0
    index_j = 1
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
    variable = nl_strain_zz
    index_i = 2
    index_j = 2
  [../]

  [./stress_xx]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_xx
    index_i = 0
    index_j = 0
  [../]
  [./stress_xy]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_xy
    index_i = 0
    index_j = 1
  [../]
  [./stress_yy]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_yy
    index_i = 1
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
  [./strain]
    type = ComputePlaneSmallStrain
  [../]
  [./stress]
    type = ComputeLinearElasticStress
  [../]
[]

[Postprocessors]
  [./disp_mag_diff]
    type = DisplacementL2Difference
    function_0 = '7.280419e-4*(x+0.5)/((x+0.5)^2+(y+0.5)^2)*(13.735473*((x+0.5)^2+(y+0.5)^2)+6.153846+(2-0.448026/((x+0.5)^2+(y+0.5)^2))*(4*(x+0.5)*(x+0.5)/((x+0.5)^2+(y+0.5)^2)-3))'
    function_1 = '7.280419e-4*(y+0.5)/((x+0.5)^2+(y+0.5)^2)*(-4.120642*((x+0.5)^2+(y+0.5)^2)-2.153846+(2-0.448026/((x+0.5)^2+(y+0.5)^2))*(3-4*(y+0.5)*(y+0.5)/((x+0.5)^2+(y+0.5)^2)))'
  [../]
  [./disp_mag_exact]
    type = FunctionsL2Norm
    function_0 = '7.280419e-4*(x+0.5)/((x+0.5)^2+(y+0.5)^2)*(13.735473*((x+0.5)^2+(y+0.5)^2)+6.153846+(2-0.448026/((x+0.5)^2+(y+0.5)^2))*(4*(x+0.5)*(x+0.5)/((x+0.5)^2+(y+0.5)^2)-3))'
    function_1 = '7.280419e-4*(y+0.5)/((x+0.5)^2+(y+0.5)^2)*(-4.120642*((x+0.5)^2+(y+0.5)^2)-2.153846+(2-0.448026/((x+0.5)^2+(y+0.5)^2))*(3-4*(y+0.5)*(y+0.5)/((x+0.5)^2+(y+0.5)^2)))'
  [../]
  [./disp_mag_error]
    type = RatioPostprocessor
    dividend = disp_mag_diff
    divisor = disp_mag_exact
  [../]
  [./stress_xx_diff]
    type = RankTwoComponentL2Difference
    rank_two_tensor = stress
    #function = '10000*(1-0.224013/((x+0.5)^2+(y+0.5)^2)*(1.5-3*(y+0.5)^2/((x+0.5)^2+(y+0.5)^2)-(0.336019/((x+0.5)^2+(y+0.5)^2)-1)*(8*(y+0.5)^4/((x+0.5)^2+(y+0.5)^2)^2-8*(y+0.5)^2/((x+0.5)^2+(y+0.5)^2)+1)))'
    function = '10000*(1-0.224013/((x+0.5)^2+(y+0.5)^2)*(1.5*cos(2*atan2(y+0.5,x+0.5))+(1-0.336019/((x+0.5)^2+(y+0.5)^2))*cos(4*atan2(y+0.5,x+0.5))))'
    index_i = 0
    index_j = 0
  [../]
  [./stress_xx_exact]
    type = FunctionsL2Norm
    #function_0 = '10000*(1-0.224013/((x+0.5)^2+(y+0.5)^2)*(1.5-3*(y+0.5)^2/((x+0.5)^2+(y+0.5)^2)-(0.336019/((x+0.5)^2+(y+0.5)^2)-1)*(8*(y+0.5)^4/((x+0.5)^2+(y+0.5)^2)^2-8*(y+0.5)^2/((x+0.5)^2+(y+0.5)^2)+1)))'
    function_0 = '10000*(1-0.224013/((x+0.5)^2+(y+0.5)^2)*(1.5*cos(2*atan2(y+0.5,x+0.5))+(1-0.336019/((x+0.5)^2+(y+0.5)^2))*cos(4*atan2(y+0.5,x+0.5))))'
  [../]
  [./stress_xx_error]
    type = RatioPostprocessor
    dividend = stress_xx_diff
    divisor = stress_xx_exact
  [../]
  [./stress_yy_diff]
    type = RankTwoComponentL2Difference
    rank_two_tensor = stress
    #function = '-10000*0.224013/((x+0.5)^2+(y+0.5)^2)*(0.5-(y+0.5)^2/((x+0.5)^2+(y+0.5)^2)+(0.336019/((x+0.5)^2+(y+0.5)^2)-1)*(8*(y+0.5)^4/((x+0.5)^2+(y+0.5)^2)^2-8*(y+0.5)^2/((x+0.5)^2+(y+0.5)^2)+1))'
    function = '10000*(-0.224013/((x+0.5)^2+(y+0.5)^2)*(0.5*cos(2*atan2(y+0.5,x+0.5))-(1-0.336019/((x+0.5)^2+(y+0.5)^2))*cos(4*atan2(y+0.5,x+0.5))))'
    index_i = 1
    index_j = 1
  [../]
  [./stress_yy_exact]
    type = FunctionsL2Norm
    #function_0 = '-10000*0.224013/((x+0.5)^2+(y+0.5)^2)*(0.5-(y+0.5)^2/((x+0.5)^2+(y+0.5)^2)+(0.336019/((x+0.5)^2+(y+0.5)^2)-1)*(8*(y+0.5)^4/((x+0.5)^2+(y+0.5)^2)^2-8*(y+0.5)^2/((x+0.5)^2+(y+0.5)^2)+1))'
    function_0 = '10000*(-0.224013/((x+0.5)^2+(y+0.5)^2)*(0.5*cos(2*atan2(y+0.5,x+0.5))-(1-0.336019/((x+0.5)^2+(y+0.5)^2))*cos(4*atan2(y+0.5,x+0.5))))'
  [../]
  [./stress_yy_error]
    type = RatioPostprocessor
    dividend = stress_yy_diff
    divisor = stress_yy_exact
  [../]
  [./stress_xy_diff]
    type = RankTwoComponentL2Difference
    rank_two_tensor = stress
    #function = '-10000*0.224013/((x+0.5)^2+(y+0.5)^2)*((x+0.5)*(y+0.5)/((x+0.5)^2+(y+0.5)^2)+(1-0.336019/((x+0.5)^2+(y+0.5)^2))*4*((x+0.5)^2-(y+0.5)^2)*(x+0.5)*(y+0.5)/((x+0.5)^2+(y+0.5)^2)^2)'
    function = '10000*(-0.224013/((x+0.5)^2+(y+0.5)^2)*(0.5*sin(2*atan2(y+0.5,x+0.5))+(1-0.336019/((x+0.5)^2+(y+0.5)^2))*sin(4*atan2(y+0.5,x+0.5))))'
    index_i = 0
    index_j = 1
  [../]
  [./stress_xy_exact]
    type = FunctionsL2Norm
    #function_0 = '-10000*0.224013/((x+0.5)^2+(y+0.5)^2)*((x+0.5)*(y+0.5)/((x+0.5)^2+(y+0.5)^2)+(1-0.336019/((x+0.5)^2+(y+0.5)^2))*4*((x+0.5)^2-(y+0.5)^2)*(x+0.5)*(y+0.5)/((x+0.5)^2+(y+0.5)^2)^2)'
    function_0 = '10000*(-0.224013/((x+0.5)^2+(y+0.5)^2)*(0.5*sin(2*atan2(y+0.5,x+0.5))+(1-0.336019/((x+0.5)^2+(y+0.5)^2))*sin(4*atan2(y+0.5,x+0.5))))'
  [../]
  [./stress_xy_error]
    type = RatioPostprocessor
    dividend = stress_xy_diff
    divisor = stress_xy_exact
  [../]
[]

[Executioner]
  type = Transient
  solve_type = PJFNK
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-11

  petsc_options = '-snes_ksp_ew'
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package -snes_lag_jacobian'
  petsc_options_value = 'lu superlu_dist 10'

  end_time = 1
  automatic_scaling = true
[]

[Outputs]
  exodus = true
[]
