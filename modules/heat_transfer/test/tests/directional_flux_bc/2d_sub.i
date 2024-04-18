[Mesh]
  [planet]
    type = ConcentricCircleMeshGenerator
    has_outer_square = false
    radii = 1
    num_sectors = 10
    rings = 2
    preserve_volumes = false
  []
  [subdomain]
    type = SubdomainBoundingBoxGenerator
    input = planet
    bottom_left = '0 -1 0'
    top_right = '1 1 0'
    block_id = 2
  []

  [moon]
    type = ConcentricCircleMeshGenerator
    has_outer_square = false
    radii = 0.5
    num_sectors = 8
    rings = 2
    preserve_volumes = false
  []
  [combine]
    type = CombinerGenerator
    inputs = 'subdomain moon'
    positions = '0 0 0 -1.5 -0.5 0'
  []
[]

[GlobalParams]
  illumination_flux = '1 1 0'
[]

[Variables]
  [u]
  []
  [v]
  []
[]

[Kernels]
  [diff_u]
    type = Diffusion
    variable = u
  []
  [dt_u]
    type = TimeDerivative
    variable = u
  []

  [diff_v]
    type = Diffusion
    variable = v
  []
  [dt_v]
    type = TimeDerivative
    variable = v
  []
[]

[BCs]
  [flux_u]
    type = DirectionalFluxBC
    variable = u
    boundary = outer
  []

  [flux_v]
    type = DirectionalFluxBC
    variable = v
    boundary = outer
    self_shadow_uo = shadow
  []
[]

[Postprocessors]
  [ave_v_all]
    type = SideAverageValue
    variable = v
    boundary = outer
  []
  [ave_v_exposed]
    type = ExposedSideAverageValue
    variable = v
    boundary = outer
    block = 2
    self_shadow_uo = shadow
  []
[]

[UserObjects]
  [shadow]
    type = SelfShadowSideUserObject
    boundary = outer
    execute_on = INITIAL
  []
[]

[Executioner]
  type = Transient
  dt = 0.1
  num_steps = 2
[]

[Outputs]
  exodus = true
[]
