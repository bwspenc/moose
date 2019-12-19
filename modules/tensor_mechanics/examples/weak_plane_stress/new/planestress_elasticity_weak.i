[GlobalParams]
  displacements = 'disp_x disp_y'
  out_of_plane_strain = nl_strain_zz
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
    function = '7.280419e-4*(x+0.5)/((x+0.5)^2+(y+0.5)^2)*(13.735473*((x+0.5)^2+(y+0.5)^2)+6.153846+(2-0.448026/((x+0.5)^2+(y+0.5)^2))*(4*(x+0.5)*(x+0.5)/((x+0.5)^2+(y+0.5)^2)-3))'
  [../]
  [./right_y]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = 2
    function = '7.280419e-4*(y+0.5)/((x+0.5)^2+(y+0.5)^2)*(-4.120642*((x+0.5)^2+(y+0.5)^2)-2.153846+(2-0.448026/((x+0.5)^2+(y+0.5)^2))*(3-4*(y+0.5)*(y+0.5)/((x+0.5)^2+(y+0.5)^2)))'
  [../]
  [./top_x]
    type = FunctionDirichletBC
    variable = disp_x
    boundary = 3
    function = '7.280419e-4*(x+0.5)/((x+0.5)^2+(y+0.5)^2)*(13.735473*((x+0.5)^2+(y+0.5)^2)+6.153846+(2-0.448026/((x+0.5)^2+(y+0.5)^2))*(4*(x+0.5)*(x+0.5)/((x+0.5)^2+(y+0.5)^2)-3))'
  [../]
  [./top_y]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = 3
    function = '7.280419e-4*(y+0.5)/((x+0.5)^2+(y+0.5)^2)*(-4.120642*((x+0.5)^2+(y+0.5)^2)-2.153846+(2-0.448026/((x+0.5)^2+(y+0.5)^2))*(3-4*(y+0.5)*(y+0.5)/((x+0.5)^2+(y+0.5)^2)))'
  [../]
[]

[Modules/TensorMechanics/Master]
  [plane_stress]
    planar_formulation = WEAK_PLANE_STRESS
    strain = SMALL
    generate_output = 'stress_xx stress_yy stress_zz strain_xx strain_yy strain_zz vonmises_stress'
    extra_vector_tags = 'ref'
  []
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

  end_time = 1.0
  automatic_scaling = true
[]

[Outputs]
  exodus = true
  csv = true
[]
