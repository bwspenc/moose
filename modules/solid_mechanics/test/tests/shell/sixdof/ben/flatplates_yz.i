[Mesh]
  [mesh]
    type = FileMeshGenerator
    file = flatplates_yz.e
  []
[]


[Variables]
  [disp_x]
  []
  [disp_y]
  []
  [disp_z]
  []
  [rot_x]
  []
  [rot_y]
  []
  [rot_z]
  []
[]

[AuxVariables]
  [react_disp_x]
  []
  [react_disp_y]
  []
  [react_disp_z]
  []
  [react_rot_x]
  []
  [react_rot_y]
  []
  [react_rot_z]
  []
[]

[BCs]
  [yz_fix_y]
    type = DirichletBC
    variable = disp_y
    boundary = '6'
    value = 0.0
  []
  [yz_fix_z]
    type = DirichletBC
    variable = disp_z
    boundary = '1'
    value = 0.0
  []
  [yz_fix_x]
    type = DirichletBC
    variable = disp_x
    boundary = '6'
    value = 0.0
  []
  [yz_fix_rot_z]
    type = DirichletBC
    variable = rot_z
    boundary = '6'
    value = 0.0
  []
  [yz_pull_y]
    type = DirichletBC
    variable = disp_y
    boundary = '8'
    value = 0.01
  []
[]

#[DiracKernels]
#  [point1]
#    type = ConstantPointSource
#    variable = disp_x
#    point = '1 0 1'
#    value = -2.5 # P = 10
#  []
#[]

#[NodalKernels]
#  [fx]
#    type = UserForcingFunctionNodalKernel
#    boundary = '3 4'
#    function = 1
#    variable = 'disp_x'
#  []
#[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  line_search = 'none'
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-8
  dt = 1.0
  dtmin = 1.0
  end_time = 1.0
[]

[Kernels]
  [solid_disp_x]
    type = ADStressDivergenceShell2
    component = 0
    variable = disp_x
    save_in = react_disp_x
    through_thickness_order = SECOND
  []
  [solid_disp_y]
    type = ADStressDivergenceShell2
    component = 1
    variable = disp_y
    save_in = react_disp_y
    through_thickness_order = SECOND
  []
  [solid_disp_z]
    type = ADStressDivergenceShell2
    component = 2
    variable = disp_z
    save_in = react_disp_z
    through_thickness_order = SECOND
  []
  [solid_rot_x]
    type = ADStressDivergenceShell2
    component = 3
    variable = rot_x
    save_in = react_rot_x
    through_thickness_order = SECOND
    penalty = 1e-6
  []
  [solid_rot_y]
    type = ADStressDivergenceShell2
    component = 4
    variable = rot_y
    save_in = react_rot_y
    through_thickness_order = SECOND
    penalty = 1e-6
  []
  [solid_rot_z]
    type = ADStressDivergenceShell2
    component = 5
    variable = rot_z
    save_in = react_rot_z
    through_thickness_order = SECOND
    penalty = 1e-6
  []
[]

[Materials]
  [elasticity]
    type = ADComputeIsotropicElasticityTensorShell
    youngs_modulus = 1e6
    poissons_ratio = 0.0
    through_thickness_order = SECOND
  []
  [strain]
    type = ADComputeIncrementalShellStrain2
    displacements = 'disp_x disp_y disp_z'
    rotations = 'rot_x rot_y rot_z'
    thickness = 0.01
    through_thickness_order = SECOND
  []
  [stress]
    type = ADComputeShellStress2
    through_thickness_order = SECOND
  []
  [stress_yy_0]
    type = RankTwoCartesianComponent
    rank_two_tensor = global_stress_t_points_0
    property_name = 'stress_yy_0t'
    index_i = 1
    index_j = 1
    outputs = all
  []
  [stress_yy_1]
    type = RankTwoCartesianComponent
    rank_two_tensor = global_stress_t_points_1
    property_name = 'stress_yy_1t'
    index_i = 1
    index_j = 1
    outputs = all
  []
[]

[Postprocessors]
  [ydisp1]
    type = PointValue
    point = '0 1 0'
    variable = disp_y
  []
  [ydisp2]
    type = PointValue
    point = '0 1 1'
    variable = disp_y
  []
  [yreact_right]
    type = NodalSum
    boundary = 8
    variable = react_disp_y
  []
  [yreact_left]
    type = NodalSum
    boundary = 6
    variable = react_disp_y
  []
[]

[Outputs]
  exodus = true
[]
