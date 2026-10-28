# Dimensionless environment verification ONLY: not a calibrated Ni-Al model.
# Checks conserved phase-field dynamics, automatic free-energy derivatives,
# periodic boundaries, PETSc nonlinear solves, MPI and Exodus output.
[Mesh]
  [domain]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 24
    ny = 24
    xmax = 24
    ymax = 24
  []
[]

[Variables]
  [c]
    [InitialCondition]
      type = SmoothCircleIC
      x1 = 12
      y1 = 12
      radius = 5
      int_width = 3
      invalue = 0.9
      outvalue = 0.1
    []
  []
  [mu]
  []
[]

[Kernels]
  [chemical_potential]
    type = SplitCHParsed
    variable = c
    w = mu
    f_name = F
    kappa_name = kappa_c
  []
  [diffusion]
    type = SplitCHWRes
    variable = mu
    mob_name = M
  []
  [conservation]
    type = CoupledTimeDerivative
    variable = mu
    v = c
  []
[]

[BCs]
  [Periodic]
    [xy]
      auto_direction = 'x y'
    []
  []
[]

[Materials]
  [constants]
    type = GenericConstantMaterial
    prop_names = 'M kappa_c'
    prop_values = '1 1'
  []
  [free_energy]
    type = DerivativeParsedMaterial
    property_name = F
    coupled_variables = c
    expression = '4*(c-0.1)^2*(c-0.9)^2'
    derivative_order = 2
  []
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  scheme = implicit-euler
  dt = 0.05
  num_steps = 10
  nl_rel_tol = 1e-9
  nl_abs_tol = 1e-11
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_type'
  petsc_options_value = 'lu mumps'
[]

[Postprocessors]
  [mean_c]
    type = ElementAverageValue
    variable = c
    execute_on = 'initial timestep_end'
  []
[]

[Outputs]
  file_base = results/environment_smoke/ch
  exodus = true
  csv = true
  execute_on = 'initial timestep_end'
[]
