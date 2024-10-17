[Mesh]
  type = GeneratedMesh
  dim=2
  nx=20
  ny=20
[]

[Problem]
  solve = false
[]

[AuxVariables]
  [av1]
  []
  [rad]
  []
[]

[Functions]
  [rad]
    type = DistanceFromPoints
    point_coordinates = '0.5 0.5 0
                         0.8 0.8 0'
  []
  [xy]
    type = PiecewiseLinear
    xy_data = '0.0  0
               0.1  1
               0.2  0
               10.0 0'
  []
  [af]
    type = AdapterFunction
    t = rad
    function = xy
  []
[]

[AuxKernels]
  [av1]
    type = FunctionAux
    variable = av1
    function = af
  []
  [rad]
    type = FunctionAux
    variable = rad
    function = rad
  []
[]

[Executioner]
  type = Steady
[]

[Outputs]
  exodus = true
[]
