# FGMFoam

![OpenFOAM v2506](https://img.shields.io/badge/OpenFOAM-v2506-blue?logo=data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAyNCAyNCI+PHBhdGggZmlsbD0id2hpdGUiIGQ9Ik0xMiAyQzYuNDggMiAyIDYuNDggMiAxMnM0LjQ4IDEwIDEwIDEwIDEwLTQuNDggMTAtMTBTMTcuNTIgMiAxMiAyeiIvPjwvc3ZnPg==)
![License: GPLv3](https://img.shields.io/badge/License-GPLv3-green.svg)
![Institution](https://img.shields.io/badge/Institution-TU%2Fe-red)

## Overview

FGMFoam is an OpenFOAM combustion solver using Flamelet Generated Manifolds (FGM) [1]. It is developed by Stijn Schepers at Eindhoven University of Technology and built on top of OpenFOAM-com (ESI/OpenCFD branch).

FGM is a tabulated chemistry approach that reduces the computational cost of combustion simulations by pre-computing the flame chemistry on a low-dimensional manifold, parameterised by a small set of control variables such as a progress variable and enthalpy. At runtime, all thermochemical quantities are retrieved from the pre-computed table rather than solved through detailed chemistry, making industrial-scale combustion simulations feasible.

FGMFoam is designed as a modular framework: all FGM-related functionality — field management, control variable transport equations, and table lookup — is self-contained in a single library with no complex cross-dependencies. Multiple FGM models can coexist in the same solver and are selected at runtime through the input files, making it straightforward to extend the framework with new FGM models for your own research without modifying the rest of the codebase.

## Prerequisites

| Requirement | Version / Notes |
|---|---|
| OpenFOAM-com (ESI/OpenCFD) | **v2506** (developed and tested on this version; later versions may work) |
| MPI | Any standard implementation (OpenMPI, MPICH) — required for parallel runs and shared-memory table loading |
| C++ standard | C++14 or later (as required by OpenFOAM-com) |
| CMake / wmake | Provided by OpenFOAM — no separate installation needed |

> **Note:** This solver is built on the **OpenFOAM-com** (ESI/OpenCFD) branch, not the OpenFOAM.org (Foundation) branch. The two branches are not interchangeable.

## Build

OpenFOAM must be sourced before building. Then source the local `bashrc` to set `LIB_LOCAL_SRC`, and run `Allwmake`:

```bash
source /path/to/openfoam/etc/bashrc   # source your OpenFOAM installation
source solver/bashrc                  # sets LIB_LOCAL_SRC=${PWD}/src
cd solver && ./Allwmake
```

`Allwmake` runs `wclean` then `wmake` in this order (order matters due to dependencies):
1. `src/lookUp/` → `$(FOAM_USER_LIBBIN)/liblookUp`
2. `src/combustionModel/` → `$(FOAM_USER_LIBBIN)/libcombustionModelFGM`
3. `src/lookUpBoundaryCondition/` → `$(FOAM_USER_LIBBIN)/liblookUpBoundaryCondition`
4. `applications/FGMFoam/` → `FGMFoam` executable

To rebuild a single library without a full clean, run `wmake` directly in the relevant subdirectory (e.g., `wmake libso src/combustionModel/`).

## Running test cases

Each test case has a `run` script that manages the full workflow:

```bash
cd testcases/laminarFlamelet && ./run
cd testcases/laminarSlotBurner && ./run
```

The `run` scripts:
1. Clean any prior results (`./clean`)
2. Copy `../../tables/database.fgm` into `constant/lookUp/`
3. Substitute placeholders in `0.orig/` and `system.orig/` via `sed`, then copy to `0/` and `system/`
4. Run `blockMesh`, `setFields`, `decomposePar`, `mpirun FGMFoam -parallel`, `reconstructPar`

Logs are written to `log.blockMesh`, `log.setFields`, `log.decomposePar`, `log.FGMFoam`, `log.reconstructPar`.

> **Note:** The number of processors is controlled by the `nProcs` variable at the top of the `run` script, and is set to 4 by default.

> **Note:** The included test cases use a methane–air FGM table at an equivalence ratio of φ = 0.7. To simulate other fuels or conditions, replace `tables/database.fgm` with a table generated for your flame configuration.

## Architecture

### Layer overview

```
FGMFoam (solver app)
    └── CombustionModel<psiReactionThermo>        [OpenFOAM interface]
            └── FGM<ReactionThermo>               [selector — reads combustionProperties]
                    └── FGM_PM_UL_NA              [concrete model: premixed + non-adiabatic effects]
                            └── FGMTable          [MPI shared-memory table manager]
                                    └── FGMlib.c  [pure-C reader + 2D interpolation]
```

### `src/lookUp` — table I/O and caching

- **`FGM/FGMlib.c`+`.h`**: Pure C library. `readFGM()` parses a `.fgm` ascii file into a `FGM` struct (pressure, control-variable grid sizes, variable names, raw float data). `lookupFGM_2D()` performs bilinear interpolation over two control variables.
- **`Table/FGMTable.C`+`.H`**: C++ singleton factory. `FGMTable::getFGMInstance(fileName)` loads the `.fgm` table once per node using an MPI shared-memory window (`MPI_Win_allocate_shared`), serialises/deserialises the `FGM` struct across ranks so only rank 0 reads the file. Returns a shared `FGM*` pointer usable by all local MPI ranks.

### `src/combustionModel` — FGM model hierarchy

- **`combustionModel/`** and **`CombustionModel/`**: Copied from OpenFOAM with no modifications; provide the base class and template selector infrastructure.
- **`FGM/FGM/FGM.C`+`.H`**: Selector. Reads `constant/combustionProperties` (`FGMCoeffs { flameType; transportModel; adiabatic; subGridModel; }`) and delegates to the matching concrete model via OpenFOAM's run-time selection.
- **`FGM/FGM_PM_UL_NA/`**: Currently the only concrete model — premixed flame with non-adiabatic effects and unity Lewis transport. Control variables are `Yc` (progress variable) and `ht` (total enthalpy). Each time step `correct()` calls `solve()` (transports Yc and ht) then `update()` (looks up T, mu, cp, lambda, rho, species mass fractions, HeatRelease from the FGM table and writes them back to the thermo fields).

### `src/lookUpBoundaryCondition` — custom boundary condition

`lookUpEnthalpyFvPatchScalarField` (`type lookUpEnthalpy`): fixed-value BC for `ht` that evaluates enthalpy at the patch by looking up T and cp from the FGM table.

### `applications/FGMFoam` — solver

PIMPLE-based compressible reactive solver. Per time step:
1. **Momentum equation** (`UEqn.H`)
2. **Combustion** — `combustion->correct()` (Yc and ht transport + FGM table lookup)
3. **Pressure corrector loop** (`pEqn.H` or `pcEqn.H` for consistent PIMPLE)
4. **Turbulence correction** — turbulent transport of the control variables is included via the gradient diffusion hypothesis, but no sub-grid scale turbulence-chemistry interaction (TCI) model is included yet.

### FGM table (`tables/database.fgm`)

ASCII table read by `FGMlib`. The format uses keyword-delimited sections:

```text
[PRESSURE]     
101325        ← operating pressure in Pa
[DIMENSION]    
2             ← number of control variables
[DATASIZE]     
250 100       ← grid points per control variable (25000 rows total)
[VARIABLES]    
16            ← number of variables, followed by their names
[END]
[DATA]       
              ← one row per table point, each row containing all Nvar values
```

Each row in `[DATA]` contains 16 space-separated floats corresponding to the variables in order: `CV1`, `CV2`, `SOURCE_CV1`, `OH`, `H2`, `H2O`, `CH4`, `CO2`, `CO`, `O2`, `DENSITY`, `TEMPERATURE`, `CP`, `CONDUCTIVITY`, `VISCOSITY`, `HEATRELEASE`. The table is stored in row-major order with CV1 (progress variable Yc) as the fast index and CV2 (total enthalpy ht) as the slow index.

The included `database.fgm` is a 250 × 100 table for a methane–air flame at φ = 0.7 and p = 101325 Pa, generated with Chem1D. Any one-dimensional flame solver can be used to generate compatible tables, provided the output is converted to the format above.

## Adding a new FGM model

1. Copy `solver/src/combustionModel/FGM/FGM_PM_UL_NA/` to `solver/src/combustionModel/FGM/YOURMODEL/`
2. Rename all files and every occurrence of `FGM_PM_UL_NA` to `YOURMODEL`
3. Add `#include "YOURMODEL.H"` in `solver/src/combustionModel/FGM/FGM/FGM.H`
4. Add `FGM/YOURMODEL/YOURMODELs.C` to `solver/src/combustionModel/Make/files`
5. Rebuild `src/combustionModel/`

## Known limitations

- **Transport models:** Only premixed flames with unity Lewis transport (`FGM_PM_UL_NA`) are currently implemented. More advanced transport models accounting for preferential and differential diffusion effects, essential for accurate hydrogen combustion modelling, are under active development.
- **Turbulence-chemistry interaction:** Turbulent simulations are supported — turbulent transport of the control variables is included via the gradient diffusion hypothesis. However, no sub-grid scale (SGS) turbulence-chemistry interaction (TCI) model is included, meaning the FGM table is accessed with instantaneous filtered control variables. Turbulent simulations are therefore limited to well-resolved LES or DNS where subgrid fluctuations are small. A SGS TCI model is planned for a future release.
- **Dimensionality:** The FGM table lookup is currently limited to two control variables (2D interpolation). Extension to higher-dimensional tables is not yet supported.
- **MPI exit crash (cosmetic):** When running in parallel, the solver may print `malloc_consolidate(): unaligned fastbin chunk detected` and MPI abort messages at the end of a completed simulation. This is a known issue originating in OpenFOAM's static destructor ordering during MPI shutdown and **does not affect simulation results or output**. The job will report `Primary job terminated normally` before the error messages appear.
- **OpenFOAM branch:** Only the OpenFOAM-com (ESI/OpenCFD) branch is supported. Compatibility with the OpenFOAM.org (Foundation) branch is not guaranteed.


## Citation

If you use FGMFoam in your research, please cite the following paper, in which the solver was applied to thermo-diffusively unstable lean premixed hydrogen–air flames:

```bibtex
@article{Schepers2025,
  author  = {S. N. J. Schepers and J. A. van Oijen},
  title   = {{FGM} modeling of thermo-diffusive unstable lean premixed hydrogen--air flames},
  journal = {Combustion and Flame},
  volume  = {280},
  year    = {2025},
  month   = oct,
  doi     = {10.1016/j.combustflame.2025.113983}
}
```

> **Note:** Hydrogen combustion capabilities described in the paper are not yet included in this public release and will be added in a future update.

## License

This project is licensed under the GNU General Public License v3.0 — see the [LICENSE](https://github.com/SNJSchepers/FGMFoam/blob/main/LICENSE) file for details.

## References

[1] J. A. van Oijen, A. Donini, R. J. Bastiaans, J. H. ten Thije Boonkkamp, L. P. de Goey, "State-of-the-art in premixed combustion modeling using flamelet generated manifolds," *Progress in Energy and Combustion Science*, 57, 30–74, 2016. https://doi.org/10.1016/j.pecs.2016.07.001
