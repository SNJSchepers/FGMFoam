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
/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.
	
\*---------------------------------------------------------------------------*/

#ifndef FGMTable_H
#define FGMTable_H

#include "FGMlib.h"
#include <mpi.h>
#include <cstring>
#include <fstream>
#include <iostream>
#include "IOstreams.H"
#include "error.H"
#include <unordered_map>
#include <unordered_set>
#include <string>


namespace Foam
{
    class FGMTable
    {
    private:

        static std::unordered_map<std::string, FGM*>*    fgmInstances_;
        static std::unordered_map<std::string, MPI_Win>* shmWins_;
        static std::unordered_map<std::string, char*>*   shmBuffers_;
        static std::unordered_map<std::string, int>*     shmBufferSizes_;
        static std::unordered_set<FGM*>*                 ownedInstances_;

        // Private constructor to prevent instantiation
        FGMTable();

        // Serialization helpers
        static int serializeFGM(const FGM* fgm, char** buffer, int* bufferSize);
        static FGM* deserializeFGM(const char* buffer, int bufferSize);

    public:
        FGMTable(const FGMTable&) = delete;
        void operator=(const FGMTable&) = delete;

        // Main API
        static FGM* getFGMInstance(const std::string& fileName);
        static void cleanup();

    };
}

#endif // FGMTable_H

 
