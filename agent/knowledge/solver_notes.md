# FGMFoam usage notes

Curated, authoritative facts about the FGMFoam solver and the laminar slot
burner case that are not obvious from the source alone. These files are indexed
into the agent's knowledge base — add or correct notes here to teach the agent.

## FGMFoam solution algorithm (what it solves, and in what order)

FGMFoam is a compressible PIMPLE (merged PISO–SIMPLE) solver, built on
rhoPimpleFoam. Per time step its loop
(`solver/applications/FGMFoam/FGMFoam.C:145-209`) does, in order:

1. **Continuity for the density**: on the first PIMPLE iteration the standard
   continuity equation `ddt(rho) + div(phi) = 0` is solved for `rho` from the
   current mass fluxes (`#include "rhoEqn.H"`, `FGMFoam.C:183-186`; the file is
   OpenFOAM's own, from `$FOAM_SRC/finiteVolume/cfdTools/compressible/`, not part
   of this repository). The momentum and scalar equations below use this density,
   which keeps their transport consistent with the fluxes.
2. **Momentum predictor**: the momentum equation is solved for the velocity with
   the current pressure gradient (`UEqn.H`, `solve(UEqn == -fvc::grad(p))`).
3. **`combustion->correct()`**: the FGM model, which is
   `FGM_PM_UL_NA::correct()` = `solve()` then `update()`
   (`solver/src/combustionModel/FGM/FGM_PM_UL_NA/FGM_PM_UL_NA.C:481-484`):
   - `solve()` solves the **two control-variable transport equations**: the
     progress variable `Yc` (convection–diffusion with a tabulated reaction
     source `sourceYc`) and the enthalpy `ht` (convection–diffusion, no source)
     (`FGM_PM_UL_NA.C:312-357`).
   - `update()` loops over the cells and looks up the FGM table at the control
     variables (Yc, ht) via `lookupFGM_2D(...)` to set the temperature, the
     compressibility `psi` (= table density / table pressure), transport
     properties (`mu`, `alpha`), species mass fractions, heat release, and the
     `Yc` source term (`FGM_PM_UL_NA.C:360-479`).
4. **Pressure corrector loop**: the pressure equation is solved and the fluxes
   and velocity are corrected (`pEqn.H`, or `pcEqn.H` for the consistent variant).
   This pressure equation is the compressible continuity equation written for
   `p` (with `rho = psi*p`). After it, `rhoEqn.H` is solved once more so that
   `compressibleContinuityErrs.H` can report the continuity error, and the density
   is then set to the thermodynamic value `rho = thermo.rho()` = `psi*p`, with
   `psi` from the FGM table (`pEqn.H:4,131`; `FGM_PM_UL_NA.C:410`).
5. Turbulence correction (`turbulence->correct()`), a no-op for the laminar case.

The governing equations are therefore: **continuity** (solved in two forms: as
the pressure equation for `p`, and as `rhoEqn.H` for `rho`), **momentum** (`U`),
and the two **control-variable transport equations** (`Yc` and `ht`). The `ht`
equation is the energy equation, in total-enthalpy form; the temperature is not
solved for but looked up from the table at (Yc, ht), and species are not
transported either: they come from the table too. The final density of each time step comes from the table and
the pressure (`rho = psi*p`); the solved continuity density keeps the transport
of U, Yc and ht conservative within the step.

The thermochemical state — density, temperature, transport, species — is NOT
computed from an ideal-gas equation of state or finite-rate chemistry; it is read
from the FGM table as a function of the two control variables (Yc, ht). Do not
describe FGMFoam as a generic combustion solver that solves species/EOS equations.

## FGM modelling assumptions (what FGM does and does NOT assume)

FGM (Flamelet-Generated Manifold) is a TABULATED-CHEMISTRY method. The assumption
is that the thermochemical state of the flame lies on a low-dimensional MANIFOLD
parameterised by a few control variables — here two: the progress variable `Yc`
and the enthalpy `ht`. The manifold is precomputed offline from one-dimensional
laminar premixed flamelet solutions with detailed chemistry and transport; at run
time the solver only transports the control variables and reads everything else
(temperature, density, transport properties, species, heat release, reaction
source) from the table. The core assumption is the FLAMELET assumption: the local
flame structure resembles that of a 1-D laminar flame, so the full
high-dimensional composition space collapses onto the manifold.

FGM does **NOT** assume chemical or "thermodynamic" equilibrium. The progress
variable parameterises FINITE-RATE chemical progress from unburnt (`Yc = 0`) to
fully burnt; the tabulated reaction source `sourceYc` is exactly the finite-rate
chemistry captured from the flamelets. Do NOT describe FGM as assuming local
thermodynamic equilibrium (LTE) or equilibrium chemistry — that is wrong; the
whole point of the progress variable is to track non-equilibrium flame progress.

The shipped model is the `FGM_PM_UL_NA` variant, selected by the four `FGMCoeffs`
switches (`flameType premixed; transportModel unityLewis; adiabatic no;
subGridModel none;`). Its specific assumptions are therefore:

- **Premixed** combustion (`flameType premixed`) — fuel and oxidiser are mixed
  before burning. The table is built at a single equivalence ratio (phi = 0.7);
  mixture fraction is NOT a control variable.
- **Unity Lewis number** (`transportModel unityLewis`, Le = 1) — all species and
  heat are assumed to diffuse at the same rate (no preferential/differential
  diffusion); the laminar diffusivity of `Yc` and `ht` is the thermal diffusivity
  `alpha`.
- **Non-adiabatic** (`adiabatic no`) — heat loss/gain is allowed, which is why
  enthalpy `ht` is carried as a SECOND control variable (an adiabatic model would
  fix enthalpy and need only `Yc`). This captures the cooled/heated walls of the
  slot burner.
- **No sub-grid model** (`subGridModel none`) — no SGS turbulence–chemistry
  closure (consistent with the laminar/DNS/well-resolved-LES use described above).

## Diffusion model: unity Lewis number (Le = 1), no differential diffusion

The diffusion (transport) model of FGMFoam's only FGM model, `FGM_PM_UL_NA`, is
**unity Lewis number**: `transportModel unityLewis` in `FGMCoeffs`, the "UL" in
the model's name. All species and heat diffuse at the same rate, so there is no
preferential or differential diffusion and no Soret effect. In practice this means
the control variables `Yc` and `ht` both diffuse with the **thermal diffusivity**
`alpha = lambda/cp`, which the table update sets in every cell
(`alphaCells[celli] = variables_[ilambda]/variables_[icp]`,
`solver/src/combustionModel/FGM/FGM_PM_UL_NA/FGM_PM_UL_NA.C:408`).

The diffusion terms in the two transport equations
(`FGM_PM_UL_NA.C:336-356`) are

    - fvm::laplacian(alpha + this->turbulence().mut()/Sct_, Yc_)
    - fvm::laplacian(alpha + this->turbulence().mut()/Sct_, ht_)

so the effective diffusivity is the laminar `alpha` (unity Lewis) plus, only when
a turbulence model is active, a turbulent part `mut/Sct` (gradient-diffusion
hypothesis, turbulent Schmidt number `Sct` default 0.7). For a laminar case
`mut` is zero and the diffusion is purely unity-Lewis molecular diffusion.

Models with non-unity Lewis numbers (preferential and differential diffusion,
needed for hydrogen flames) are not in this release; a new transport model would
be a new FGM model class selected through `transportModel`.

## The FGM update step: what the table lookup writes into the thermo fields

`FGM_PM_UL_NA::update()` (`solver/src/combustionModel/FGM/FGM_PM_UL_NA/FGM_PM_UL_NA.C:360-479`)
runs after the Yc/ht transport, every time step. It loops over **every cell, and
then every boundary face**, looks up the table at that point's (Yc, ht) with
`lookupFGM_2D`, and writes the result into fields it finds by name:

- `T` (temperature), `mu` (viscosity), `alpha` = lambda/cp (thermal
  diffusivity, which under unity Lewis is also the diffusivity of Yc and ht);
- `thermo:psi` = table density / table pressure, so that `rho = psi*p`;
- the species mass fractions `CH4`, `H2`, `H2O`, `OH`, `O2`, `CO2`, `CO`, and
  `HeatRelease`;
- `sourceYc`, the source term of the Yc equation.

Three rules are applied on the way: a **negative `sourceYc` is clipped to 0**
(cells and faces); `sourceYc` is **set to 0 on `wall` patches** (no reaction at a
wall); and on a boundary the **temperature is only overwritten where `T` does not
fix its value** (`!pT.fixesValue()`), so a prescribed wall or inlet temperature
is kept. The thermo fields are written directly (the model takes them by name and
const_casts them), which is why renaming a table column or a field breaks this.

## Boundary conditions of the laminar slot burner case

The field boundary conditions live in the case's time-zero field files
(`testcases/laminarSlotBurner/0.orig/{U,p,T,Yc,ht}`), NOT in
`constant/polyMesh/boundary` — that file only defines the mesh patch geometry
(patch names and their geometric type: `inlet`/`outlet` are `patch`, `plate` is a
`wall`, `symm_top`/`symm_bottom` are `symmetry`, `front`/`back` are `empty`). To
answer "what BC is used for X", read the `0.orig/<field>` file, not the polyMesh
boundary file. The case is LAMINAR, so there are no turbulence fields and no wall
functions. Build_case fills the numeric values (the ALLCAPS placeholders below).

The seven patches are `inlet`, `outlet`, `plate`, `symm_top`, `symm_bottom`,
`front`, `back`. On every field, `symm_top`/`symm_bottom` are `symmetry` and
`front`/`back` are `empty` (the case is 2D). The remaining (inlet/outlet/plate)
boundary conditions per field:

- **U** (`0.orig/U`): inlet `fixedValue` (uniform `(INLETVELOCITY 0 0)`, axial
  inflow); outlet `zeroGradient`; plate `noSlip`.
- **p** (`0.orig/p`): inlet `zeroGradient`; outlet `waveTransmissive` (a
  non-reflecting outflow, with `fieldInf 101325`, `gamma 1.25`, `lInf 0.01`);
  plate `zeroGradient`. (A plain `fixedValue 101325` outlet is present but
  commented out.)
- **T** (`0.orig/T`): inlet `fixedValue` (`INLETTEMPERATURE`); outlet
  `zeroGradient`; plate `fixedValue` (`PLATETEMPERATURE`, the hot wall).
- **Yc** (`0.orig/Yc`, the progress variable): inlet `fixedValue` (the unburnt
  value, 0 for the cold premixed inflow); outlet `zeroGradient`; plate
  `zeroGradient`.
- **ht** (`0.orig/ht`, the enthalpy control variable): inlet and plate use the
  custom `lookUpEnthalpy` BC (see the next section — it sets the enthalpy from the
  prescribed wall/inlet temperature via the FGM table); outlet `zeroGradient`.

So temperature is imposed directly (fixedValue) on the inlet and the hot plate,
and the matching enthalpy at those patches is derived from it by `lookUpEnthalpy`;
the outlet is a non-reflecting pressure outflow with zeroGradient scalars.

## The laminarFlamelet case: 1D freely propagating flame

`testcases/laminarFlamelet/` is the second example case, and it is physically a
different problem from the slot burner, not a variant of it.

- **Geometry:** a 1D domain (patches `open`, `closed`, `sides`), with `sides`
  of type `empty` — the usual OpenFOAM way of making a mesh one- or
  two-dimensional.
- **No inlet and no walls in the flow sense:** `U` is `inletOutlet` on `open`
  and `noSlip` on `closed`; `T`, `Yc` and `ht` are **`zeroGradient` on every
  patch**. There is no `fixedValue` inlet, no plate temperature, and — unlike
  the slot burner — **no `lookUpEnthalpy` boundary condition** anywhere.
- **The mixture starts quiescent:** `internalField` for `U` is `(0 0 0)`. The
  flame is started by `setFields`, which sets `Yc` to the burnt value
  (`FINALPROGRESSVARIABLE`) inside a `boxToCell` region from x = 1.5 cm to
  5 cm, leaving the rest at the unburnt value. The flame front then propagates
  freely into the fresh mixture.
- **Only conditions are parameterised:** `INITIALTEMPERATURE` (note: *not*
  `INLETTEMPERATURE`, which is the slot burner's spelling), `INITIALPRESSURE`,
  `INITIALPROGRESSVARIABLE`, `FINALPROGRESSVARIABLE`, `INITIALENTHALPY`, plus
  the run controls `SIMULATIONTIME`, `SIMULATIONDELTATIME`,
  `SIMULATIONWRITETIME` and `NPROCS`.
- **Typical time step is smaller than the slot burner's:** the shipped `run`
  script uses `deltaT = 2e-6 s` (vs 1e-5 s for the slot burner), because the
  propagating front must stay resolved in time.

## How long to run the slot burner: flow-through time

A steady slot-burner flame is only a result once the fresh gas has swept the
domain and flushed out the initial condition. The flow-through time is the axial
domain length divided by the inlet velocity: for the shipped geometry, 28 mm at
0.2 m/s gives about 0.14 s. The pre-run check warns when endTime is shorter than
one flow-through, and notes when it is fewer than about 3, after which a steady
flame is normally converged. This is conservative — burnt gas downstream moves
faster than the inlet gas. It does not apply to the laminarFlamelet, which has no
inflow; there the flame propagates into gas at rest.

## How the time step is chosen

The solver runs at a FIXED step (`adjustTimeStep no`; the `maxCo` entry in
controlDict is inert), so a bad choice cannot be corrected during the run. The
agent computes one from the mesh and the flow whenever the user does not give
one:

    dt = CFL * dx / U,  with CFL = 0.1

The velocity is **not** the inlet velocity. Gas leaves the flame about
rho_u/rho_b times faster than it enters (~6x for methane-air at phi = 0.7), and
that burnt-gas velocity sets the Courant number:

- **laminarSlotBurner**: U = inlet_velocity * rho_u/rho_b, dx = flame_thickness/10
  (the real cell size). Costs nothing extra.
- **laminarFlamelet**: nothing but the flame moves the gas, so U = S_L *
  rho_u/rho_b from a 1D free-flame solve (unity Lewis, no Soret, matching the FGM model), dx = 10 um (fixed by its blockMeshDict).

This was checked against real runs: the formula gives 5.07e-6 s for the slot
burner at the shipped conditions, and a run at 5e-6 s held a settled Courant of
0.11 against the 0.1 target.

## The Courant number spikes at ignition — judge the settled value

Both cases are ignited as a step function in Yc, which is a discontinuity rather
than a flame. Relaxing it spikes the Courant number for a handful of steps. A
slot burner run that sat at 0.11 for 6000 steps reported a maximum of 4.11,
entirely from step 3. `analyse_run` reports both the overall maximum and the
maximum after the first 20 steps, and only a **settled** value above 1 means the
step was genuinely too large.

## Recommended time step, and why a hotter mixture needs a smaller one

Each case ships with run-control values in its own `run` script, and those are
what the agent uses when the user does not give them:

- `laminarSlotBurner`: `sim_time` 0.1 s, `sim_timestep` **1e-5 s**
- `laminarFlamelet`: `sim_time` 0.1 s, `sim_timestep` **2e-6 s** (five times
  smaller — the front propagates freely and must stay resolved in time)

Both were chosen at an unburnt temperature of **300 K**.

The step is **fixed**: `controlDict` sets `adjustTimeStep no`, and the `maxCo`
entry there is inert (it only matters for adaptive stepping). Nothing reduces
the step at run time, so a step that is too large for the conditions does not
self-correct — the run diverges instead.

A higher unburnt temperature therefore needs a **smaller** step. The laminar
flame speed of methane-air rises steeply with unburnt temperature (roughly
T^1.5-2) and the flame front gets thinner, so the same step resolves less of the
front. As a conservative starting point the agent suggests scaling the step by
(300 K / T)^2 — e.g. 2e-6 s becomes about 5e-7 s at 600 K — and recommends
checking the Courant number printed in `log.FGMFoam`, reducing further if it
approaches 1. This is guidance, not a derived stability criterion.

## How the flame is measured after a run

The agent reports these from the written fields; both use a temperature
threshold of **900 K** by default to mark the flame front. Methane's
autoignition temperature is about 810 K (537 C), so 900 K sits just above it,
well inside the reaction zone and far from both the unburnt (300 K) and burnt
(~1840 K at phi = 0.7) states. The exact value matters less than using the same
one consistently.

**laminarFlamelet — laminar flame speed.** The front position is tracked over
the written times and fitted against time. Two corrections matter:

- The speed is taken **relative to the unburnt gas**, not in the lab frame. The
  burnt gas expands and pushes the fresh mixture towards the open end, so the
  lab-frame front speed is larger than the flame speed by that gas velocity.
- The case is ignited by `setFields` as a **step function** in Yc, which is a
  discontinuity and not a flame structure. The first part of the run is the
  profile relaxing into a real flame, and its apparent speed during that time is
  meaningless, so it is excluded before fitting.

**laminarSlotBurner — flame length.** The flame anchors on the plate and closes
on the **slot centreline (y = h)** downstream of it. The length is the axial
distance from the plate's **trailing edge** to the first point on the slot
centreline reaching the threshold. The plate's position comes from the other
symmetry line (y = 0, the wake behind the plate), where it leaves a gap in the
sampled points. The wake is NOT where the flame length is measured: it fills with
burnt gas just behind the plate at any inlet velocity. Three outcomes replace the
number rather than qualify it:

- **flashback** — the gas is above the threshold upstream of the plate, i.e. the
  flame has propagated into the inlet channel;
- **blow-off** — nothing downstream of the holder reaches the threshold, so no
  flame is anchored;
- **longer_than_domain** — the wake is burnt but the slot centreline is cold up
  to the outlet, or the tip is within 5% of it: the flame is anchored, but the
  domain length, not the flame, sets the answer. It is NOT blow-off. The fix is a
  larger `outlet_length`. At the default geometry the length grows about 4.5 mm
  per 0.1 m/s of inlet velocity and the tip reaches the default outlet at about
  0.45 m/s.

## The slot burner's geometry is parameterised

`testcases/laminarSlotBurner/system.orig/blockMeshDict` carries four
placeholders, so the burner shape is a case parameter and not a dictionary edit:

- `PLATETHICKNESS` -> `l_plate`, the streamwise thickness of the plate that
  holds the flame;
- `HOLEPITCH` -> `h = pitch/2`; the domain is **one half-slot**, symmetric about
  y = 0 (the wake centreline behind the plate) and y = h (the slot centreline);
- `HOLESIZE` -> `h_plate = h - size/2`, so this is the open height of the slot;
- `OUTLETLENGTH` -> `l_out`, the domain length downstream of the plate
  (`outlet_length`); the flame tip must fall inside it;
- `FLAMETHICKNESS` -> `lF`, which sets the cell size `dx = lF/10` **and nothing
  else** — resolution and shape are independent knobs.

The defaults reproduce the original hardcoded mesh exactly (`l_plate = 5 lF`,
`h = 10 lF`, `h_plate = 4 lF`, `l_out = 30 lF`), giving the same 43000 cells, so
parameterising it changed no result. The bash `run` script sets the same values
at the top, so the case is still usable without the agent.

The **laminarFlamelet's mesh is not parameterised** — its geometry cannot be
varied or swept.

Because the flame-length measurement locates the plate from the gap in the
sampled wake-line points, and the slot centreline from the mesh's upper edge,
rather than from fixed numbers, it follows a changed geometry automatically.

## Which cases the agent can set up

The agent builds either template via `build_case`'s `case_type`:
`laminarSlotBurner` (the default) or `laminarFlamelet`. It varies conditions
and, for the slot burner only, the five geometry parameters above; the
boundary-condition *types* and the numerical schemes always come from the
template unchanged. To change boundary-condition types or the flamelet's mesh,
edit the template case by hand; the agent will not do it.

The flamelet needs only `T_inlet` (the unburnt-gas temperature, also accepted as
`T_initial`/`T_unburnt`/`T_0`), `sim_time` and `sim_timestep`. Asking it for an
inlet velocity or a plate temperature is a category error: that case has
neither.

## Enthalpy boundary condition: lookUpEnthalpy (temperature -> enthalpy inversion)

The transported energy variable in FGMFoam is the enthalpy control variable `ht`,
but at walls and inlets the user knows a TEMPERATURE, not an enthalpy. The custom
boundary condition `lookUpEnthalpy` (class `lookUpEnthalpyFvPatchScalarField`,
`solver/src/lookUpBoundaryCondition/lookUp/lookUpEnthalpyFvPatchScalarField.C`)
bridges that gap: it is a **fixedValue** patch field for `ht` whose value is the
enthalpy that makes the FGM table return the prescribed patch temperature. In the
shipped laminar slot burner it is used on the `inlet` and `plate` patches of the
`ht` field (`0.orig/ht`); the `h0` entry is the initial enthalpy guess.

How `updateCoeffs()` sets the value (lines 304-378), per patch face:

1. Read the target temperature `T0` from the patch `T` field, and set the two
   control variables from the patch — the progress variable `Yc` and the current
   enthalpy.
2. Look up the FGM table (`lookupFGM_2D`) at those control variables to get the
   table temperature `T1` and heat capacity `cp`. This lookup is plain bilinear
   interpolation between the stored grid points of the precomputed
   `database.fgm` table: numbers computed beforehand from 1D flame solutions.
3. **Iterate to converge the temperature** (NOT to retry a failed lookup): nudge
   the enthalpy control variable by `dha = 0.1*(T0 - T1)*cp` and re-look up the
   table, repeating until `|T0 - T1| <= 1 K` (the residual `Tresidual = 1.0`), to
   a maximum of 1000 iterations (it warns if it hits the cap).
4. Fix the patch enthalpy to the converged value.

So it is a relaxed fixed-point inversion of the tabulated relation `T(ht)`: given
the boundary temperature it finds the matching enthalpy. The factor `cp` and the
`0.1` under-relaxation make it a damped Newton step. Do NOT describe it as looking
up T and cp to "compute" the enthalpy directly, and do not describe the loop as
retrying failed table lookups — it converges the enthalpy so the table temperature
equals the prescribed boundary temperature.

## Numerical setup: schemes, linear solvers, time stepping (laminar slot burner)

Grounded in the case's system files (`testcases/laminarSlotBurner/system.orig/`).

Discretisation (`fvSchemes`):

- **Time** (`ddtSchemes` default): `Euler` — first-order implicit Euler (transient),
  chosen for stability. The time term is NOT "Gauss linear" — Gauss linear is a
  spatial scheme; the ddt scheme here is Euler.
- **Convection** (`divSchemes`): `default none` (so every convection term must be
  named explicitly). `div(phi,U)` is `Gauss linear` (second-order central);
  `div(phi,Yc)` is `Gauss limitedLinear 0` and **`ht` uses the same scheme as
  `Yc`** — a limited-linear (bounded) scheme for the transported scalars. The
  viscous-stress term `div((rho*nuEff)*dev2(T(grad(U))))` is `Gauss linear`.
- **Gradients** `Gauss linear`; **diffusion** (`laplacianSchemes`) `Gauss linear
  orthogonal`; **interpolation** `linear`; **surface-normal gradient**
  (`snGradSchemes`) `orthogonal`.

Linear solvers (`fvSolution`, `solvers`):

- `p` and `rho`: **`PCG`** (preconditioned conjugate gradient) with **`DIC`**
  preconditioner, tol 1e-8 / relTol 1e-3; the `(rho|p)Final` corrector tightens to
  tol 1e-9 / relTol 0.
- `U`, `Yc`, `ht`: **`PBiCGStab`** with **`DILU`** preconditioner, tol 1e-8 /
  relTol 1e-3; the `Final` variant tol 1e-9 / relTol 0.
- Pressure–velocity coupling (`PIMPLE`): `momentumPredictor yes`, `nCorrectors 1`,
  `nNonOrthogonalCorrectors 0`, `transonic no`. `nOuterCorrectors` is unset
  (defaults to 1), so PIMPLE runs effectively in PISO mode (one outer corrector
  per step).

Time stepping (`controlDict`):

- The time step is **FIXED**: `adjustTimeStep no`, so `deltaT` is held CONSTANT for
  the whole run. Its value is set by build_case's `sim_timestep` (and editable
  later with change_settings). The `maxCo 1.0` and `maxDeltaT 1` entries are
  present but INACTIVE because adjust-time-step is off — do NOT describe the run as
  Courant-limited or adaptive.
- Output: `writeControl adjustableRunTime`, `writeInterval` = end time /
  sim_write_steps, `runTimeModifiable true`.

## Progress variable definition: Yc = Y_CO2 / M_CO2

The progress variable Yc is the specific mole number of CO2:

    Yc = Y_CO2 / M_CO2

where Y_CO2 is the CO2 mass fraction and M_CO2 its molar mass (≈ 44.0 kg/kmol),
so Yc is in kmol/kg. This is the definition used to generate the FGM table and to
set the case's initial/boundary values (it is computed this way in
`agent/cantera_tool.py`).

It is based on CO2 ONLY (one species divided by its molar mass) and is NOT
normalised to 0–1. For the shipped methane/air table at phi = 0.7, Yc = 0 in the
unburnt mixture (no CO2) and Yc ≈ 2.44e-3 kmol/kg at chemical equilibrium (burnt).
Do NOT define Yc as a normalised (0–1) sum of CO / CO2 / H2 / H2O mass fractions —
that generic FGM progress-variable definition is not the one used here.

## Turbulence: well-resolved LES/DNS work; no SGS-TCI closure; RANS not accurate

FGMFoam is NOT restricted to laminar flow. It creates a standard OpenFOAM
compressible turbulence model (`solver/applications/FGMFoam/createFields.H:78-91`),
and turbulent transport of the control variables is included through a
gradient-diffusion term — the turbulent diffusivity `turbulence().mut()/Sct_` is
added to the diffusion of the Yc and ht equations (`FGM_PM_UL_NA.C:341,353`). What
is NOT implemented is a sub-grid-scale turbulence–chemistry interaction (SGS-TCI)
model: `createFields.H:89-91` prints "WARNING: No Sub-grid scale model
incorporated in this solver" whenever the turbulence model is not laminar
(`Stokes`).

Consequences (so a flat "the solver can only do laminar" is WRONG):

- **DNS and well-resolved (fine) LES: supported and accurate** — when the flame is
  resolved on the grid, the unresolved sub-grid fluctuations and the associated
  turbulence–chemistry interaction are negligible, so no SGS-TCI closure is needed.
- **Coarse LES and RANS: not accurate** — they rely on the unresolved fluctuations
  that this solver does not model. RANS in particular needs a model for the
  fluctuations (TCI), which is not implemented.
- **Laminar: supported.**

Separately, the two SHIPPED example cases (`laminarSlotBurner`, `laminarFlamelet`)
ARE laminar — their `constant/turbulenceProperties` set `simulationType laminar` —
but that is a property of those cases, not a limitation of the solver. Distinguish
"what the solver can do" from "how the example cases are configured".

## The FGM table is always constant/lookUp/database.fgm

The FGM lookup table is always read from `constant/lookUp/database.fgm`. There is
no configurable table path anywhere in the case — there is no `tablePath` entry
and no file that points to the table by name. Do not claim that a table path is
configured. The path is hard-coded in the model's constructor in
`solver/src/combustionModel/FGM/FGM_PM_UL_NA/FGM_PM_UL_NA.C`:
`fgm_ = FGMTable::getFGMInstance("constant/lookUp/database.fgm");`. Using another
table means copying it to that path (the case's `run` script and the agent both
do this) or changing that line and recompiling.

## combustionProperties selects the FGM model, not a table path

`constant/combustionProperties` sets `combustionModel FGM` and configures the FGM
model through the `FGMCoeffs` dictionary. `FGMCoeffs` contains exactly four word
switches, read by the solver in `solver/src/combustionModel/FGM/FGM/FGM.C` (lines
63-66):

    FGMCoeffs
    {
        flameType      premixed;
        transportModel unityLewis;
        adiabatic      no;
        subGridModel   none;
    }

These four switches select the concrete FGM model variant — here
premixed + unityLewis + adiabatic=no + subGridModel=none picks the `FGM_PM_UL_NA`
model — and the solver errors on any unknown value. The only optional extra the
solver reads is a `turbulenceCoeffs` subdictionary holding a turbulent Schmidt
number `Sct` (default 0.7, `FGM_PM_UL_NA.C:287-291`), absent from the example case.
Those four switches (plus the optional turbulenceCoeffs) are the only keys in this
file — do not invent or attribute any others to it.

## thermophysicalProperties is an unused placeholder

`constant/thermophysicalProperties` is NOT used for the thermochemistry. Its own
header states that the models in the file are not used and that all
thermophysical parameters are overwritten by the FGM combustion model. It holds
only generic placeholder values (e.g. a single `pureMixture` with `molWeight
28.9`, `Cp 1000`) so the OpenFOAM framework initialises without error. It has no
`fuel`, `oxidant`, or `inert` entries — do not infer the fuel, oxidiser, or
mixture from it. All thermochemistry comes from the FGM table.

## Operating pressure must match the FGM table

The FGMFoam solver itself can run at any operating pressure. The constraint is the
FGM table: it tabulates the thermochemistry at a single pressure, recorded in its
header under the `[PRESSURE]` key (in Pa). The operating/initial pressure `p_0`
must equal that table pressure. The FGM table shipped with the repository
(`tables/database.fgm`, copied into each case's `constant/lookUp/database.fgm`)
was generated at `[PRESSURE] 101325` Pa (1 atm). Running at a different pressure
(for example 10 bar) is therefore physically inconsistent with this table. The
solver handles any pressure, but the tabulated chemistry would not match, so do
NOT build a case at a different pressure without a matching regenerated table:
use 101325 Pa and tell the user a different pressure requires regenerating the
FGM table at that pressure first.

## The FGM table is fixed at equivalence ratio 0.7

The current FGM table (`database.fgm`) is generated for a single equivalence
ratio of phi = 0.7. Equivalence ratio (mixture fraction) is NOT one of the
table's control variables / dimensions: the table is two-dimensional, tabulated
over the progress variable and enthalpy (CV1 and CV2 in the header) at one fixed
mixture composition. The simulation must therefore use phi = 0.7. A different
equivalence ratio is not supported by this table — the tabulated chemistry would
be inconsistent with it, so you would need to regenerate the FGM table at the
desired equivalence ratio. When setting up a case, default the equivalence ratio
to 0.7 and tell the user it is fixed by the table; do not build a case at a
different phi without a matching regenerated table.

## The FGM table header (database.fgm)

`database.fgm` is an ASCII file whose header (the lines before the `[DATA]`
block) records its metadata: `[PRESSURE]` (generation pressure, Pa),
`[DIMENSION]` (manifold dimension / number of control variables), `[DATASIZE]`
(grid resolution), and `[VARIABLES]` (the tabulated quantities — the control
variables CV1 and CV2, the source term SOURCE_CV1, species mass fractions such as
CH4, CO2, CO, O2, H2O, OH and H2, and DENSITY, TEMPERATURE, CP, CONDUCTIVITY,
VISCOSITY, HEATRELEASE). Only this header is indexed for the agent; the bulk
numeric data after `[DATA]` is not.

## Viewing results / post-processing (what is written, what the agent can do)

FGMFoam writes standard OpenFOAM TIME directories in the case directory (e.g.
`testcases/laminarSlotBurner_T300/0.001/`), reconstructed from the per-rank
`processor*/` dirs by `reconstructPar`. Each written time holds the transported
fields `U`, `p`, `T`, `Yc`, `ht`, plus the quantities the FGM model writes from the
table (all `AUTO_WRITE` in `FGM_PM_UL_NA.C`): the species mass fractions `CH4`,
`O2`, `CO2`, `CO`, `H2`, `H2O`, `OH`, the heat-release rate `HeatRelease`, and the
`Yc` reaction source `sourceYc`. So temperature and the major species ARE written
and available to look at.

To view them, point the user to **ParaView**: `paraFoam -case <caseDir>`, or create
an empty `<caseDir>/case.foam` (or `.OpenFOAM`) file and open it in ParaView. For
sampling / slicing / derived quantities use the standard OpenFOAM `postProcess`
utility (or `sample` / function objects). There is NO special FGMFoam
post-processor — do NOT invent tool names like "OFPost" or commands like
"plotAllFields"/"-suffix .p"; only the standard OpenFOAM utilities and ParaView
exist here. And note there is NO `postProcessing/` directory unless function
objects are enabled (they are not — the controlDict functions block is commented
out).

Keep this guidance BRIEF and command-level (the two options above). Do NOT write a
step-by-step ParaView GUI/menu walkthrough, do NOT tell the user to open
`log.FGMFoam` or the `0.orig` directory in ParaView (the log is text; `0.orig` is
the initial template, not results), and do NOT name ParaView "widgets" (e.g. a
"Temperature widget") for fields — those are fabrications. ParaView loads the whole
case and exposes every written field for colouring/slicing; just point the user to
`paraFoam`/the `.foam` file.

The AGENT ITSELF cannot read, plot, analyse, or interpret the numeric field data
(there is no tool for it); it can only report run status and the log tail via
run_status. For looking at fields, always direct the user to ParaView or the
standard OpenFOAM utilities — never claim to plot or analyse results yourself.

The table is a precomputed, structured grid of values read from this ASCII file
and queried by direct bilinear interpolation in `lookupFGM_2D`
(`solver/src/lookUp/FGM/FGMlib.c`). Every value the solver gets is interpolated
from stored numbers that were computed beforehand from 1D flamelet solutions;
"Flamelet-Generated Manifold" is the name of that tabulated chemistry. Describe
the lookup as "interpolation in a precomputed table".

## Table lookup: bilinear interpolation, and what happens outside the table range

Every value the solver takes from the FGM table (temperature, density, viscosity,
conductivity, cp, species, heat release, the Yc source term) comes from one C
function: `lookupFGM_2D` in `solver/src/lookUp/FGM/FGMlib.c`. It interpolates
**bilinearly** in the two control variables:

- CV2 (`ht`) is normalised over its range, which is the same for every row;
- CV1 (`Yc`) is normalised over the range of *that* `ht` row, which is itself
  interpolated, because the reachable Yc range depends on the enthalpy;
- the four surrounding table points are weighted (`w1m`, `w1p`, `w2m`, `w2p`).

**Outside the table range it extrapolates linearly; it does not clip.** The cell
index is clamped to the first or last cell (`fgm_max(0, i)`, `fgm_min(N-1, i)`),
but the weights are not, so a point past the edge gets a weight above 1 or below 0.
The result is a linear extrapolation from the edge cell. The function always
returns `EXIT_SUCCESS`: there is no error, warning or abort for an out-of-range
lookup. Values far outside the table can therefore be unphysical (for example a
negative source term, which `FGM_PM_UL_NA` then clips to 0). This is why the
agent's pre-run check refuses initial `Yc`/`ht` values outside the table's ranges.

## Table columns are matched by name — renaming a column breaks the lookup

The solver finds each quantity in the FGM table by the **name** of its column in
the `[VARIABLES]` header, not by position. The constructor in
`solver/src/combustionModel/FGM/FGM_PM_UL_NA/FGM_PM_UL_NA.C` loops over the
table's variable names and compares each one with `std::strcmp` against fixed
strings: `SOURCE_CV1`, `TEMPERATURE`, `CP`, `CONDUCTIVITY`, `VISCOSITY`, `DENSITY`,
`HEATRELEASE`, and the species `CH4`, `H2`, `H2O`, `OH`, `O2`, `CO2`, `CO`, storing
the matching column index (`isourceYc`, `iT`, ...). These strings are the contract
with the table format. A column that is renamed, for example `TEMPERATURE` to `T`,
is simply never matched: its index stays at its initial value of -1, no error is
raised at construction, and the lookups that use it read the wrong memory. Keep
the names exactly as above when generating a new table.
