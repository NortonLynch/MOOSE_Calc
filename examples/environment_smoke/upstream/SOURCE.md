# Official MOOSE KKS input

`kks_example_split.i` is copied unchanged from the INL MOOSE conda package
`moose=2026.10.14=mpich`, source commit `c140bfd5a6`.

Upstream path: `modules/phase_field/test/tests/KKS_system/kks_example_split.i`.

Source: https://github.com/idaholab/moose/blob/c140bfd5a6/modules/phase_field/test/tests/KKS_system/kks_example_split.i

MOOSE licensing: LGPL 2.1; https://github.com/idaholab/moose/blob/c140bfd5a6/LICENSE
Copyright information: https://github.com/idaholab/moose/blob/c140bfd5a6/COPYRIGHT

The adjacent `../kks_smoke.i` derived copy adds a concentration-average postprocessor,
project-local outputs and uses PETSc ASM/LU with the correct `-pc_type` spelling.
It retains the required subdomain factor shift and removes only the unused
global factor-shift option.
The original input is a generic toy problem, not Ni-Al thermodynamics.
