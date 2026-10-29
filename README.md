# FGMFoam

![OpenFOAM v2606](https://img.shields.io/badge/OpenFOAM-v2506%20%7C%20v2606-blue?logo=data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAyNCAyNCI+PHBhdGggZmlsbD0id2hpdGUiIGQ9Ik0xMiAyQzYuNDggMiAyIDYuNDggMiAxMnM0LjQ4IDEwIDEwIDEwIDEwLTQuNDggMTAtMTBTMTcuNTIgMiAxMiAyeiIvPjwvc3ZnPg==)
![License: GPLv3](https://img.shields.io/badge/License-GPLv3-green.svg)
![Institution](https://img.shields.io/badge/Institution-TU%2Fe-red)

## Overview

FGMFoam is an OpenFOAM combustion solver using Flamelet Generated Manifolds (FGM) [1]. It is developed by Stijn Schepers at Eindhoven University of Technology and built on top of OpenFOAM-com (ESI/OpenCFD branch).

FGM is a tabulated chemistry approach that reduces the computational cost of combustion simulations by pre-computing the flame chemistry on a low-dimensional manifold, parameterised by a small set of control variables such as a progress variable and enthalpy. At runtime, all thermochemical quantities are retrieved from the pre-computed table rather than solved through detailed chemistry, making industrial-scale combustion simulations feasible.

FGMFoam is designed as a modular framework: all FGM-related functionality (field management, control variable transport equations, and table lookup) is self-contained in a single library with no complex cross-dependencies. Multiple FGM models can coexist in the same solver and are selected at runtime through the input files, making it straightforward to extend the framework with new FGM models for your own research without modifying the rest of the codebase.

An optional **conversational agent** (`agent/`) can answer questions about the
solver from its own source code and set up, check, run and analyse cases from a
plain-language description; see [The FGMFoam agent](#the-fgmfoam-agent-optional)
below and its full documentation in [`agent/README.md`](agent/README.md). The
solver does not need it.

## Contents

- [Prerequisites](#prerequisites) · [Build](#build) · [Running test cases](#running-test-cases)
- [Architecture](#architecture) · [Adding a new FGM model](#adding-a-new-fgm-model) · [Known limitations](#known-limitations)
- [The FGMFoam agent (optional)](#the-fgmfoam-agent-optional), full docs in [`agent/README.md`](agent/README.md)
- [Citation](#citation) · [License](#license) · [References](#references)

## Prerequisites

| Requirement | Version / Notes |
|---|---|
| OpenFOAM-com (ESI/OpenCFD) | **v2506** and **v2606**: developed on v2506, and both test cases build and run on v2606. Other versions may work |
| MPI | Any standard implementation (OpenMPI, MPICH), required for parallel runs and shared-memory table loading |
| C++ standard | C++14 or later (as required by OpenFOAM-com) |
| CMake / wmake | Provided by OpenFOAM, no separate installation needed |
| Python | 3.10+, **only** for the optional [agent](#the-fgmfoam-agent-optional); the solver needs no Python |

> **Note:** This solver is built on the **OpenFOAM-com** (ESI/OpenCFD) branch, not the OpenFOAM.org (Foundation) branch. The two branches are not interchangeable.

## Build

OpenFOAM must be sourced before building. Then, from inside `solver/`, source the local `bashrc` to set `LIB_LOCAL_SRC`, and run `Allwmake`:

```bash
source /path/to/openfoam/etc/bashrc   # source your OpenFOAM installation
cd solver
source bashrc                         # sets LIB_LOCAL_SRC=${PWD}/src, so run it from solver/
./Allwmake
```

`Allwmake` always exits with status 0, even when a step fails: check that `FGMFoam` exists afterwards (`which FGMFoam`), or look for `Error` in its output.

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

> **Note:** The slot burner's mesh is parameterised: `plate_thickness`, `hole_pitch`, `hole_size`, `outlet_length` and the flame thickness `lF` (which sets the cell size) are set at the top of its `run` script and substituted into `blockMeshDict`. The defaults reproduce the original geometry, so results are unchanged; varying them changes the burner shape without editing the dictionary.

> **Tip:** To vary the conditions you edit the parameters at the top of the `run` script. The optional [agent](#the-fgmfoam-agent-optional) does that substitution for you from a plain-language description, writing a fresh case alongside the template rather than editing it.

> **Note:** The included test cases use a methane–air FGM table at an equivalence ratio of φ = 0.7. To simulate other fuels or conditions, replace `tables/database.fgm` with a table generated for your flame configuration.

## Architecture

### Layer overview

```
FGMFoam (solver app)
    └── CombustionModel<psiReactionThermo>        [OpenFOAM interface]
            └── FGM<ReactionThermo>               [selector: reads combustionProperties]
                    └── FGM_PM_UL_NA              [concrete model: premixed + non-adiabatic effects]
                            └── FGMTable          [MPI shared-memory table manager]
                                    └── FGMlib.c  [pure-C reader + 2D interpolation]
```

### `src/lookUp`: table I/O and caching

- **`FGM/FGMlib.c`+`.h`**: Pure C library. `readFGM()` parses a `.fgm` ascii file into a `FGM` struct (pressure, control-variable grid sizes, variable names, raw float data). `lookupFGM_2D()` performs bilinear interpolation over two control variables.
- **`Table/FGMTable.C`+`.H`**: C++ singleton factory. `FGMTable::getFGMInstance(fileName)` loads the `.fgm` table once per node using an MPI shared-memory window (`MPI_Win_allocate_shared`), serialises/deserialises the `FGM` struct across ranks so only rank 0 reads the file. Returns a shared `FGM*` pointer usable by all local MPI ranks.

### `src/combustionModel`: FGM model hierarchy

- **Base classes**: `combustionModel`, `CombustionModel<ReactionThermo>` and `ThermoCombustion` come from OpenFOAM's own `libcombustionModels`, which the library links against. FGMFoam registers its model into OpenFOAM's run-time selection table rather than compiling its own copy of those classes. Older versions did compile such copies, so their static objects existed twice and were destroyed twice at exit, which made parallel runs end with `malloc_consolidate(): unaligned fastbin chunk detected` (exit code 134).
- **`FGM/FGM/FGM.C`+`.H`**: Selector. Reads `constant/combustionProperties` (`FGMCoeffs { flameType; transportModel; adiabatic; subGridModel; }`) and delegates to the matching concrete model via OpenFOAM's run-time selection.
- **`FGM/FGM_PM_UL_NA/`**: Currently the only concrete model: premixed flame with non-adiabatic effects and unity Lewis transport. Control variables are `Yc` (progress variable) and `ht` (total enthalpy). Each time step `correct()` calls `solve()` (transports Yc and ht) then `update()` (looks up T, mu, cp, lambda, rho, species mass fractions, HeatRelease from the FGM table and writes them back to the thermo fields).

### `src/lookUpBoundaryCondition`: custom boundary condition

`lookUpEnthalpyFvPatchScalarField` (`type lookUpEnthalpy`): fixed-value BC for `ht` that evaluates enthalpy at the patch by looking up T and cp from the FGM table.

### `applications/FGMFoam`: solver

PIMPLE-based compressible reactive solver. Per time step:
1. **Momentum equation** (`UEqn.H`)
2. **Combustion**: `combustion->correct()` (Yc and ht transport + FGM table lookup)
3. **Pressure corrector loop** (`pEqn.H` or `pcEqn.H` for consistent PIMPLE)
4. **Turbulence correction**: turbulent transport of the control variables is included via the gradient diffusion hypothesis, but no sub-grid scale turbulence-chemistry interaction (TCI) model is included yet.

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
- **Turbulence-chemistry interaction:** Turbulent simulations are supported: turbulent transport of the control variables is included via the gradient diffusion hypothesis. However, no sub-grid scale (SGS) turbulence-chemistry interaction (TCI) model is included, meaning the FGM table is accessed with instantaneous filtered control variables. Turbulent simulations are therefore limited to well-resolved LES or DNS where subgrid fluctuations are small. A SGS TCI model is planned for a future release.
- **Dimensionality:** The FGM table lookup is currently limited to two control variables (2D interpolation). Extension to higher-dimensional tables is not yet supported.
- **OpenFOAM branch:** Only the OpenFOAM-com (ESI/OpenCFD) branch is supported. Compatibility with the OpenFOAM.org (Foundation) branch is not guaranteed.

## The FGMFoam agent (optional)

> **Full documentation: [`agent/README.md`](agent/README.md)**: installation,
> model presets, the tool list, adding your own case, and the measured results.
> The design rationale is in [`agent/DESIGN.md`](agent/DESIGN.md).

`agent/` contains an optional conversational assistant for this solver: a CLI that
answers questions about FGMFoam grounded in its own source code, and sets up, runs
and monitors cases from a plain-language description. Both example cases can be
built (the 2D `laminarSlotBurner` and the 1D `laminarFlamelet`), varying the
conditions, and for the slot burner the burner geometry as well (plate
thickness, hole pitch, hole size, outlet length), since its mesh is parameterised.

```
> What is the progress variable definition?
  Yc = Y_CO2 / M_CO2: the CO2 mass fraction over its molar mass (kmol/kg), not a
  normalised 0-1 quantity. (agent/knowledge/solver_notes.md)

> Set up phi 0.7, U 0.2 m/s, inlet 300 K, plate 300 K, run 0.1 s, dt 1e-5 s
  Case written to testcases/laminarSlotBurner_agent/  (ready to run)

> Now a flamelet case at 300 K, 0.1 s, dt 2e-6
  Case written to testcases/laminarFlamelet_agent/    (ready to run)

> How did that run go?
  Reached endTime, max Courant 0.42. Laminar flame speed 0.19 m/s
  (initialisation transient excluded, R2 = 0.999).
```

It replaces the manual `sed` editing in the test cases' `run` scripts with a
validated, deterministic pipeline: the language model only decides *when* to build
or run, while parameter validation, the Cantera thermochemistry, and the case
writing are plain Python and identical no matter which model is used. Every
generated case carries an `agent-provenance.json` recording the parameters, which
of them you chose, the mechanism and the FGM table's content hash; and before the
solver is launched the case is checked against that table: conditions outside
its tabulated range, or a time step too large for the mesh, are caught in
milliseconds rather than after hours of compute. Afterwards it reports what
happened: the Courant number the run actually reached, and the **laminar flame
speed** (flamelet) or **flame length**, flashback or blow-off (slot burner). Answers are
grounded in a local retrieval index over the solver source, and for "how does X
work" questions the agent opens the actual source file before answering.

It runs on **the Mistral API (the default), the Claude API, the OpenAI API, or a
fully local open-source model** (Ollama/vLLM; nothing leaves the machine, which
matters for unpublished research). Switching takes one word, at launch or
mid-conversation, and you can simply ask:

```bash
agent/.venv/bin/python agent/main.py --model claude     # or mistral (default), openai, local
> /models                          # list the presets      (/cases lists the cases)
> /model local-3b                  # switch mid-conversation
> use Claude for this one          # or just ask, same thing
```

### Quick start

Python 3.10+, run from the repository root:

```bash
python3 -m venv agent/.venv
agent/.venv/bin/pip install -r agent/requirements.txt
```

Then pick how it should run. The default is Mistral Large 4 ("Le Chonk"), which needs a key:

```bash
export MISTRAL_API_KEY=your_key_here     # or ANTHROPIC_API_KEY / OPENAI_API_KEY
agent/.venv/bin/python agent/main.py
```

Or run it entirely on your own machine, with no key and no data leaving it
(needs [Ollama](https://ollama.com); see
[`agent/README.md` §2.3](agent/README.md#23-running-fully-locally-no-data-leaves-the-machine)
for the one-time model pull):

```bash
agent/.venv/bin/python agent/main.py --model local
```

On first run the agent builds its knowledge index over the solver source; that
takes a moment and happens once.

The two example cases are demonstrators: **your own simulations can be added as
plugins**: one folder with a short declaration of the case's parameters,
placeholders and what to measure, and they get provenance, pre-run checks,
automatic time steps, campaigns and graphs without touching the agent's code.

See [`agent/README.md`](agent/README.md) for installation, configuration, model
presets, [adding your own case](agent/README.md#5-adding-your-own-case), the full tool list, and the reliability notes;
[`agent/DESIGN.md`](agent/DESIGN.md) records the design rationale.

> **Note:** The agent is a convenience layer. The solver, the test cases, and the
> `run` scripts are fully usable without it, and it is not required to build or run
> FGMFoam.

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
  doi     = {10.1016/j.combustflame.2025.114332}
}
```

> **Note:** Hydrogen combustion capabilities described in the paper are not yet included in this public release and will be added in a future update.

## License

This project is licensed under the GNU General Public License v3.0; see the [LICENSE](https://github.com/SNJSchepers/FGMFoam/blob/main/LICENSE) file for details.

## References

[1] J. A. van Oijen, A. Donini, R. J. Bastiaans, J. H. ten Thije Boonkkamp, L. P. de Goey, "State-of-the-art in premixed combustion modeling using flamelet generated manifolds," *Progress in Energy and Combustion Science*, 57, 30–74, 2016. https://doi.org/10.1016/j.pecs.2016.07.001
