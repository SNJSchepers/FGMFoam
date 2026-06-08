# FGMFoam

![OpenFOAM v2506](https://img.shields.io/badge/OpenFOAM-v2506-blue?logo=data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAyNCAyNCI+PHBhdGggZmlsbD0id2hpdGUiIGQ9Ik0xMiAyQzYuNDggMiAyIDYuNDggMiAxMnM0LjQ4IDEwIDEwIDEwIDEwLTQuNDggMTAtMTBTMTcuNTIgMiAxMiAyeiIvPjwvc3ZnPg==)
![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)
![Institution](https://img.shields.io/badge/Institution-TU%2Fe-red)

## Overview

FGMFoam is an OpenFOAM combustion solver using Flamelet Generated Manifolds (FGM). It is developed by Stijn Schepers at Eindhoven University of Technology and built on top of OpenFOAM-com (ESI/OpenCFD branch).

FGM is a tabulated chemistry approach that reduces the computational cost of combustion simulations by pre-computing the flame chemistry on a low-dimensional manifold, parameterised by a small set of control variables such as a progress variable and enthalpy. At runtime, all thermochemical quantities are retrieved from the pre-computed table rather than solved through detailed chemistry, making industrial-scale combustion simulations feasible.

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

## Architecture

### Layer overview

```
FGMFoam (solver app)
    └── CombustionModel<psiReactionThermo>   [OpenFOAM interface]
            └── FGM<ReactionThermo>          [selector — reads combustionProperties]
                    └── FGM_PM_CH4_HL        [concrete model: premixed CH4 + heat loss]
                            └── FGMTable     [MPI shared-memory table manager]
                                    └── FGMlib.c  [pure-C reader + 2D interpolation]
```

### `src/lookUp` — table I/O and caching

- **`FGM/FGMlib.c`+`.h`**: Pure C library. `readFGM()` parses a `.fgm` binary file into a `FGM` struct (pressure, control-variable grid sizes, variable names, raw float data). `lookupFGM_2D()` performs bilinear interpolation over two control variables.
- **`Table/FGMTable.C`+`.H`**: C++ singleton factory. `FGMTable::getFGMInstance(fileName)` loads the `.fgm` table once per node using an MPI shared-memory window (`MPI_Win_allocate_shared`), serialises/deserialises the `FGM` struct across ranks so only rank 0 reads the file. Returns a shared `FGM*` pointer usable by all local MPI ranks.

### `src/combustionModel` — FGM model hierarchy

- **`combustionModel/`** and **`CombustionModel/`**: Copied from OpenFOAM with no modifications; provide the base class and template selector infrastructure.
- **`FGM/FGM/FGM.C`+`.H`**: Selector. Reads `constant/combustionProperties` (`FGMCoeffs { flameType; fuel; heatloss; subGridModel; }`) and delegates to the matching concrete model via OpenFOAM's run-time selection.
- **`FGM/FGM_PM_CH4_HL/`**: Currently the only concrete model — premixed methane flame with heat loss. Control variables are `Yc` (progress variable) and `ht` (total enthalpy). Each time step `correct()` calls `solve()` (transports Yc and ht) then `update()` (looks up T, mu, cp, lambda, rho, species mass fractions, HeatRelease from the FGM table and writes them back to the thermo fields).

### `src/lookUpBoundaryCondition` — custom boundary condition

`lookUpEnthalpyFvPatchScalarField` (`type lookUpEnthalpy`): fixed-value BC for `ht` that evaluates enthalpy at the patch by looking up T and cp from the FGM table.

### `applications/FGMFoam` — solver

PIMPLE-based compressible reactive solver. Per time step:
1. **Momentum equation** (`UEqn.H`)
2. **Combustion** — `combustion->correct()` (Yc and ht transport + FGM table lookup)
3. **Pressure corrector loop** (`pEqn.H` or `pcEqn.H` for consistent PIMPLE)
4. **Turbulence correction** — note: no sub-grid scale model for the turbulence-chemistry interaction is included yet; a warning is printed if the turbulence model is not `Laminar`

### FGM table (`tables/database.fgm`)

Binary table read by `FGMlib`. The table is indexed by two control variables (Yc, ht for the premixed CH4 model) and stores all dependent thermo and species fields. The test cases copy this table to `constant/lookUp/database.fgm` at run time.

## Adding a new FGM model

1. Copy `solver/src/combustionModel/FGM/FGM_PM_CH4_HL/` to `solver/src/combustionModel/FGM/YOURMODEL/`
2. Rename all files and every occurrence of `FGM_PM_CH4_HL` to `YOURMODEL`
3. Add `#include "YOURMODEL.H"` in `solver/src/combustionModel/FGM/FGM/FGM.H`
4. Add `FGM/YOURMODEL/YOURMODELs.C` to `solver/src/combustionModel/Make/files`
5. Rebuild `src/combustionModel/`

## Known limitations

- **Fuel support:** Only premixed methane–air flames (`FGM_PM_CH4_HL`) are currently implemented. Hydrogen combustion support is under active development; note that accurate hydrogen flame modelling requires accounting for preferential diffusion effects (non-unity Lewis numbers), which are not yet included.
- **Turbulence modelling:** No sub-grid scale (SGS) turbulence-chemistry (TCI) interaction model is included. The solver will print a warning if the turbulence model is not set to `Laminar`. SGS TCI support is planned.
- **Dimensionality:** The FGM table lookup is currently limited to two control variables (2D interpolation). Extension to higher-dimensional tables is not yet supported.
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

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.
