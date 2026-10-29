# OpenFOAM / FGMFoam naming conventions (so the source reads correctly)

Authoritative meanings of the symbols and operators used in the FGMFoam source, so
explanations do not mislabel them. Grounded in the solver
(`solver/applications/FGMFoam/createFields.H`,
`solver/src/combustionModel/FGM/FGM_PM_UL_NA/FGM_PM_UL_NA.C`).

## Fields and variables

- **`U`** — velocity (m/s). **`p`** — pressure (Pa). **`T`** — temperature (K).
- **`rho`** — density (kg/m^3); set from the thermo as `rho = thermo.rho()`.
- **`phi`** — the **mass flux**, a `surfaceScalarField` with units kg/s, defined as
  `phi = linearInterpolate(rho*U) & mesh.Sf()` (`createFields.H`). It is NOT the
  velocity. In `fvm::div(phi, X)` / `mvConvection().fvmDiv(phi, X)` it carries the
  convective transport of `X`.
- **`alpha`** — the **laminar thermal diffusivity** `thermo.alpha()` (= lambda/cp,
  units kg/(m*s)), NOT a viscosity (`createFields.H`). Under the unity-Lewis
  transport model (Le = 1) it is ALSO the laminar diffusivity of the control
  variables `Yc` and `ht`.
- **`mu`** / **`nu`** — dynamic / kinematic molecular viscosity.
  **`turbulence().mut()`** (`mut`) — the turbulent (eddy) dynamic viscosity.
- **`Sct_`** — the turbulent Schmidt number (default **0.7**; from a
  `turbulenceCoeffs` subdict). The turbulent scalar diffusivity is `mut/Sct_`, so
  the **effective diffusivity** of `Yc` and `ht` in their transport equations is
  `alpha + turbulence().mut()/Sct_` = laminar (`alpha`) + turbulent (`mut/Sct`).
- **`psi`** — the compressibility `thermo:psi` = (FGM-table density)/(FGM-table
  pressure); the density relation is `rho = psi*p`. It comes from the table, not an
  ideal-gas EOS (`FGM_PM_UL_NA.C` sets `psi = variables_[irho]/fgm_->pressure`).
- **`Yc`** — the progress variable, `Yc = Y_CO2 / M_CO2` (see the progress-variable
  note); **`ht`** — the enthalpy control variable. These two are the FGM control
  variables (the table's CV1 and CV2).
- **`sourceYc_`** — the tabulated reaction source term for `Yc`, looked up per cell
  from the FGM table and clamped to be non-negative. The `ht` equation has no
  source.

## Finite-volume operators (fvm = implicit, fvc = explicit)

- **`fvm::ddt(rho, X)`** — the TRANSIENT term d(rho*X)/dt. Not convection.
- **`fvm::div(phi, X)`** or **`mvConvection().fvmDiv(phi, X)`** — CONVECTION,
  div(phi*X), with `phi` the mass flux. `mvConvection` is a *multivariate
  convection (surface-interpolation) scheme* shared by the transported scalars — it
  is NOT the FGM table lookup.
- **`fvm::laplacian(D, X)`** — DIFFUSION, div(D*grad(X)).
- A field on the right of `==` is an explicit **source**; `fvm::Sp` / `fvm::SuSp`
  are implicit sources.
- **`fvScalarMatrix`** — the assembled linear system for a scalar transport
  equation; `.relax()` under-relaxes it and `.solve()` solves it.
- **`psiReactionThermo`** — the compressible (psi-based) reaction thermo type the
  solver is built on; here its thermodynamics are overwritten by the FGM table.
