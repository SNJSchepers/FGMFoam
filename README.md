# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

FGMFoam is an OpenFOAM combustion solver using Flamelet Generated Manifolds (FGM). It is developed by Stijn Schepers at Eindhoven University of Technology and built on top of OpenFOAM-com (ESI/OpenCFD branch).

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
1. **Density equation** (`rhoEqn.H`)
2. **Momentum equation** (`UEqn.H`)
3. **Combustion** — `combustion->correct()` (Yc and ht transport + FGM table lookup)
4. **Pressure corrector loop** (`pEqn.H` or `pcEqn.H` for consistent PIMPLE)
5. **Turbulence correction** — note: no sub-grid scale model is included; a warning is printed if the turbulence model is not `Stokes`

### FGM table (`tables/database.fgm`)

Binary table read by `FGMlib`. The table is indexed by two control variables (Yc, ht for the premixed CH4 model) and stores all dependent thermo and species fields. The test cases copy this table to `constant/lookUp/database.fgm` at run time.

## Adding a new FGM model

1. Copy `solver/src/combustionModel/FGM/FGM_PM_CH4_HL/` to `solver/src/combustionModel/FGM/YOURMODEL/`
2. Rename all files and every occurrence of `FGM_PM_CH4_HL` to `YOURMODEL`
3. Add `#include "YOURMODEL.H"` in `solver/src/combustionModel/FGM/FGM/FGM.H`
4. Add `FGM/YOURMODEL/YOURMODELs.C` to `solver/src/combustionModel/Make/files`
5. Rebuild `src/combustionModel/`
