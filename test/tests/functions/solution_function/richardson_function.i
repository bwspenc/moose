[Mesh]
  #This is the input or output mesh for whichever refinement level you want to do this
  #integration on
  file = my_mesh.e 
[]

[Problem]
  solve = false
[]

[AuxVariables]
  [dummy]
    initial_condition = 0
  []
  [disp_y_1]
  []
  [disp_y_2]
  []
  [disp_y_n]
  []
[]

[AuxKernels]
  [disp_y_1]
    type = FunctionValueAux
    variable = disp_y_1
    function = function_1
  []
  [disp_y_2]
    type = FunctionValueAux
    variable = disp_y_2
    function = function_2
  []
  [disp_y_n]
    type = FunctionValueAux
    variable = disp_y_n
    function = function_n
  []
[]

[Materials]
  [error]
    type = ParsedMaterial
    property_name = 'error'
    coupled_variables = 'disp_y_1 disp_y_2 disp_y_n'
    expression = 'disp_y_n-4*disp_y_1-3*disp_y_2'
  []
[]

[Functions]
  [function_1]
    type = SolutionFunction
    solution = solution_1
    from_variable = 'disp_y'
  []
  [function_2]
    type = SolutionFunction
    solution = solution_2
    from_variable = 'disp_y'
  []
  [function_n]
    type = SolutionFunction
    solution = solution_n
    from_variable = 'disp_y'
  []
[]

[UserObjects]
  [solution_1]
    type = SolutionUserObject
    system_variables = 'disp_x disp_y disp_z'
    #This is the output mesh with the most refined solution
    mesh = solution_1_out.e
  []
  [solution_2]
    type = SolutionUserObject
    system_variables = 'disp_x disp_y disp_z'
    #This is the output mesh with the second most refined solution
    mesh = solution_2_out.e
  []
  [solution_n]
    type = SolutionUserObject
    system_variables = 'disp_x disp_y disp_z'
    #This is the output mesh with the nth refined solution that you compute the error on
    mesh = solution_n_out.e
  []
[]

[Postprocessors]
  [richardson_error]
    type = RichardsonExtraplationError
    variable = dummy
    function_1 = function_1
    function_2 = function_2
    function_n = function_n
  []
  [rel_error]
    type = ElementL2Error
    variable = disp_y_n
    function = function_1
  []
  [richardson_error_parsed]
    type = ElementIntegralMaterialProperty
    mat_prop = error
  []
[]

[Executioner]
  type = Transient
  #Use time stepping to match your solutions. They don't have to match -- it will interpolate if
  #they don't match.
  start_time = 0.0
  end_time = 1.0
  dt = 1.0
[]

[Outputs]
  file_base = out
  csv = true
[]
