# Ni–Al lookup versus quadratic free energy

The original case in `../ni_al_coarsening/` is preserved. Both new inputs have
identical initial nodal concentration/order fields, geometry, phase definitions,
mesh, interface parameters, mobilities and reference scales. They replay the
accepted time sequence of the original run through `t*=60000`.

- `lookup.i`: nested KKS with the original 502/242-point Thermo-Calc tables.
- `quadratic_control.i`: the same nested KKS machinery with the original quadratics.
- `quadratic_matched_tail.i`: fresh `SolutionIC` restart from the step-1000 saved fields using
  all lookup-accepted times, including ten extra subdivisions near dissolution.
  This quadratic tail uses the original algebraic elimination (validated against
  the complete nested control), avoiding extra local-Newton failures at dissolution.
- `*_hermite_coefficients.csv`: interval coefficients for energy in units of
  1000 J/mol, as `a3*(c-left)^3 + a2*(c-left)^2 + a1*(c-left) + a0`.
- `baseline_time_sequence.csv`: prescribed integration times, not merely output times.
- `lookup_accepted_time_sequence.csv`: the final accepted schedule for the matched comparison.

`scripts/prepare_lookup_comparison.py` constructs cubic Hermite interpolants from
the CSV columns `x_Al`, `Gm_reduced_J_mol`, `dGm_reduced_dx_J_mol`. It generates
a balanced interval-search expression for each phase. MOOSE's
`DerivativeParsedMaterial` evaluates the selected table interval and takes
consistent derivatives; this is not a global polynomial fit. Both G and G' are
continuous; G'' is piecewise continuous. Reduced energies must use the same
affine reference for both phases. The damped local KKS solve enforces the
exported composition bounds, without clipping accepted compositions.

The finite-element integration-point phase concentrations are material
properties. Exodus outputs `c_gamma` and `c_gamma_prime` as element averages;
`c`, `eta`, and `mu` are nodal variables. The global initial fields are identical
to the baseline, even though the new thermodynamics subsequently adjusts the
local phase partition. Mobilities, not a composition-independent chemical
diffusivity rederived from the new curvature, are held fixed.

Run a fresh comparison from the project root in WSL. Install the `moose`
environment as described in the root README, then create the pinned analysis
environment and download the published reference fields (one-time setup):

```bash
source ~/miniforge3-moose/etc/profile.d/conda.sh
conda create -y -n moose-post --file environment/postprocess-linux-64.lock.txt
conda activate moose-post
python scripts/download_lookup_results.py
```

If `moose-post` already exists, skip `conda create`. The download script verifies
the archive and each file against the committed SHA256 manifest, restores the
original relative paths, and refuses to overwrite different local files.
The archive includes the original `coarsening.e` needed by the analysis.
See the [published result guide](../../results/ni_al_lookup/README.md).

Then run the complete comparison:

```bash
bash scripts/run_lookup_comparison.sh
```

The launcher generates a new timestamped result directory and refuses to
overwrite it. It completes both full runs, then replays the quadratic tail on
the lookup-accepted schedule as a time-refinement sensitivity check. The final
quadratic tail itself required two additional accepted subdivisions: the
time grids are exactly identical through `t*=40000` (both rate-fit windows),
but not throughout dissolution. Late endpoint differences must therefore be
interpreted with that remaining temporal-discretization limitation. This
restart uses the step-1000 Exodus frame. Since implicit Euler needs only the
previous solution, the saved c/eta/mu fields suffice; local KKS concentrations are
recomputed. `assemble_lookup_control.py` verifies those initial fields and joins
the unchanged original prefix and newly computed tail into `quadratic_matched.e`.
No field resampling is performed. A direct checkpoint recovery was tested but
restored the old timestepper schedule. `quadratic_matched.i` and
`prepare_matched_time_control.py` are retained only as diagnostic source files;
the launcher uses the validated `quadratic_matched_tail.i` route above.
To regenerate only this completed comparison's figures/report:

```bash
source ~/miniforge3-moose/etc/profile.d/conda.sh
conda activate moose-post
python scripts/analyze_lookup_comparison.py
```

Results are in `results/ni_al_lookup/REPORT.html`, with PNG/PDF plots, the
comparison GIF, particle-history CSVs, metrics JSON and two ParaView-readable
Exodus files. The primary pair is `lookup.e` and `quadratic_matched.e`;
`quadratic_control.e` preserves the original-schedule numerical control.
Reference time remains dimensionless `t*`, not calibrated aging
time. A four-particle 2-D simulation does not determine a 3-D LSW rate constant.

Official implementation references:
- https://mooseframework.inl.gov/source/materials/KKSPhaseConcentrationMaterial.html
- https://mooseframework.inl.gov/source/materials/DerivativeParsedMaterial.html
- https://mooseframework.inl.gov/source/timesteppers/TimeSequenceStepper.html
