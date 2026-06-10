/*---------------------------------------------------------------------------*\

███████╗ ██████╗ ███╗   ███╗███████╗ ██████╗  █████╗ ███╗   ███╗
██╔════╝██╔════╝ ████╗ ████║██╔════╝██╔═══██╗██╔══██╗████╗ ████║
█████╗  ██║  ███╗██╔████╔██║█████╗  ██║   ██║███████║██╔████╔██║
██╔══╝  ██║   ██║██║╚██╔╝██║██╔══╝  ██║   ██║██╔══██║██║╚██╔╝██║
██║     ╚██████╔╝██║ ╚═╝ ██║██║     ╚██████╔╝██║  ██║██║ ╚═╝ ██║
╚═╝      ╚═════╝ ╚═╝     ╚═╝╚═╝      ╚═════╝ ╚═╝  ╚═╝╚═╝     ╚═╝

FGMFoam - Combustion solver using Flamelet Generated Manifolds 
Developed by Stijn Schepers - Eindhoven University of Technology

FGMFoam is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published
by the Free Software Foundation, either version 3 of the License,
or (at your option) any later version.

\*---------------------------------------------------------------------------*/

#ifndef FGMLIB_H
#define FGMLIB_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <sys/stat.h>

#ifdef __cplusplus
extern "C" {
#endif

#define VAR_NAME_LENGTH 50
#define MAX_LINE_LENGTH 256

#define fgm_max(a,b) \
   ({ __typeof__ (a) _a = (a); \
      __typeof__ (b) _b = (b); \
      _a > _b ? _a : _b; })

#define fgm_min(a,b) \
   ({ __typeof__ (a) _a = (a); \
      __typeof__ (b) _b = (b); \
      _a < _b ? _a : _b; })

///@brief Structure for reading the fgm data tables
typedef struct {
    float pressure;                       ///< Pressure at which the table is constructed.
    int Ncv;                              ///< Number of control variables.
    int *Ngrid;                           ///< Number of grid points per dimension.
    int Ntotal;                           ///< Total number of grid points.
    int Nvar;                             ///< Number of variables.
    char (*varname)[VAR_NAME_LENGTH];     ///< Variable names.
    float *data;                          ///< Raw table data.
} FGM;

/// @brief Reads a full FGM table from file.
/// @param filename Path to FGM file.
/// @return Pointer to allocated FGM structure.
FGM *readFGM(const char filename[]);

/// @brief Frees memory for an FGM structure.
/// @param fgm Pointer to FGM structure.
/// @return 0 if success.
int freeFGM(FGM *fgm);

/// @brief Locates a keyword in the FGM file.
/// @param fid File pointer to the FGM file.
/// @param keyword Keyword to search for.
/// @return 0 if succes.
int locateKeyWord(FILE *fid, char keyword[]);

/// @brief Performs a 2D lookup in the FGM table.
/// @param fgm Pointer to the FGM structure
/// @param x Pointer to the control variables.
/// @param f Pointer to the output variables
int lookupFGM_2D (FGM *fgm, float *x, float *f);


#ifdef __cplusplus
}
#endif

#endif // CHEMISTRY_H
