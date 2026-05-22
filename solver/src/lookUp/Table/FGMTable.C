/*---------------------------------------------------------------------------*\

███████╗ ██████╗ ███╗   ███╗███████╗ ██████╗  █████╗ ███╗   ███╗
██╔════╝██╔════╝ ████╗ ████║██╔════╝██╔═══██╗██╔══██╗████╗ ████║
█████╗  ██║  ███╗██╔████╔██║█████╗  ██║   ██║███████║██╔████╔██║
██╔══╝  ██║   ██║██║╚██╔╝██║██╔══╝  ██║   ██║██╔══██║██║╚██╔╝██║
██║     ╚██████╔╝██║ ╚═╝ ██║██║     ╚██████╔╝██║  ██║██║ ╚═╝ ██║
╚═╝      ╚═════╝ ╚═╝     ╚═╝╚═╝      ╚═════╝ ╚═╝  ╚═╝╚═╝     ╚═╝

FGMFoam - Combustion solver using Flamelet Generated Manifolds 
Developed by Stijn Schepers - Eindhoven University of Technology

FGMFoam is developed for research and experimental purposes.
Use, modification, and redistribution are subject to the
applicable OpenFOAM license terms.

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

#include "FGMTable.H"

namespace Foam
{

    // Static member definitions
    std::unordered_map<std::string, FGM*> FGMTable::fgmInstances_;
    std::unordered_map<std::string, MPI_Win> FGMTable::shmWins_;
    std::unordered_map<std::string, char*> FGMTable::shmBuffers_;
    std::unordered_map<std::string, int> FGMTable::shmBufferSizes_;

    FGMTable::FGMTable() {}

    FGM* FGMTable::getFGMInstance(const std::string& fileName)
    {
        // Return cached instance if already loaded on this rank
        auto it = fgmInstances_.find(fileName);
        if (it != fgmInstances_.end())
        {
            return it->second;
        }

        // Check if MPI is initialized
        int initialized = 0;
        MPI_Initialized(&initialized);

        if (!initialized)
        {
            Info << "MPI not initialized. Loading table in local memory: "
                 << fileName << nl << endl;

            FGM* localFGM = readFGM(fileName.c_str());
            fgmInstances_[fileName] = localFGM;
            return localFGM;
        }

        Info << "Loading FGM table into shared memory: "
             << fileName << nl << endl;

        // Create shared-memory communicator for this node
        MPI_Comm shmComm;
        MPI_Comm_split_type
        (
            MPI_COMM_WORLD,
            MPI_COMM_TYPE_SHARED,
            0,
            MPI_INFO_NULL,
            &shmComm
        );

        int localRank = 0;
        MPI_Comm_rank(shmComm, &localRank);

        // References to per-file storage
        MPI_Win& shmWin = shmWins_[fileName];
        char*& shmBuffer = shmBuffers_[fileName];
        int& shmBufferSize = shmBufferSizes_[fileName];

        if (localRank == 0)
        {
            Pout << "Local rank 0 loading table from disk: "
                 << fileName << nl << endl;

            FGM* localFGM = readFGM(fileName.c_str());

            char* localBuffer = nullptr;
            if (serializeFGM(localFGM, &localBuffer, &shmBufferSize) != 0)
            {
                freeFGM(localFGM);
                MPI_Comm_free(&shmComm);

                FatalErrorInFunction
                    << "Failed to serialize FGM table: " << fileName
                    << exit(FatalError);
            }

            MPI_Win_allocate_shared
            (
                shmBufferSize,
                sizeof(char),
                MPI_INFO_NULL,
                shmComm,
                &shmBuffer,
                &shmWin
            );

            std::memcpy(shmBuffer, localBuffer, shmBufferSize);

            free(localBuffer);
            freeFGM(localFGM);
        }
        else
        {
            Pout << "Local rank " << localRank
                 << " accessing shared memory for table: "
                 << fileName << nl << endl;

            MPI_Aint size = 0;
            int dispUnit = 0;

            MPI_Win_allocate_shared
            (
                0,
                sizeof(char),
                MPI_INFO_NULL,
                shmComm,
                &shmBuffer,
                &shmWin
            );

            MPI_Win_shared_query
            (
                shmWin,
                0,
                &size,
                &dispUnit,
                &shmBuffer
            );

            shmBufferSize = static_cast<int>(size);
        }

        MPI_Win_fence(0, shmWin);
        MPI_Comm_free(&shmComm);

        FGM* fgm = deserializeFGM(shmBuffer, shmBufferSize);
        fgmInstances_[fileName] = fgm;

        Info << "FGM table successfully loaded into shared memory: "
             << fileName << nl << endl;

        return fgm;
    }

    int FGMTable::serializeFGM(const FGM* fgm, char** buffer, int* bufferSize)
    {
        int size =
            sizeof(float) +                              // pressure
            sizeof(int) * 3 +                            // Ncv, Ntotal, Nvar
            sizeof(int) * fgm->Ncv +                     // Ngrid
            VAR_NAME_LENGTH * fgm->Nvar +                // varname
            sizeof(float) * fgm->Ntotal * fgm->Nvar;    // data

        *bufferSize = size;
        *buffer = static_cast<char*>(malloc(size));
        if (!*buffer) return -1;

        char* ptr = *buffer;

        memcpy(ptr, &fgm->pressure, sizeof(float));
        ptr += sizeof(float);

        memcpy(ptr, &fgm->Ncv, sizeof(int));
        ptr += sizeof(int);

        memcpy(ptr, &fgm->Ntotal, sizeof(int));
        ptr += sizeof(int);

        memcpy(ptr, &fgm->Nvar, sizeof(int));
        ptr += sizeof(int);

        memcpy(ptr, fgm->Ngrid, sizeof(int) * fgm->Ncv);
        ptr += sizeof(int) * fgm->Ncv;

        for (int i = 0; i < fgm->Nvar; ++i)
        {
            memcpy(ptr, fgm->varname[i], VAR_NAME_LENGTH);
            ptr += VAR_NAME_LENGTH;
        }

        memcpy(ptr, fgm->data, sizeof(float) * fgm->Ntotal * fgm->Nvar);
        ptr += sizeof(float) * fgm->Ntotal * fgm->Nvar;

        return 0;
    }

    FGM* FGMTable::deserializeFGM(const char* buffer, int bufferSize)
    {
        if (!buffer || bufferSize <= 0) return nullptr;

        const char* ptr = buffer;
        FGM* fgm = static_cast<FGM*>(malloc(sizeof(FGM)));
        if (!fgm) return nullptr;

        memcpy(&fgm->pressure, ptr, sizeof(float));
        ptr += sizeof(float);

        memcpy(&fgm->Ncv, ptr, sizeof(int));
        ptr += sizeof(int);

        memcpy(&fgm->Ntotal, ptr, sizeof(int));
        ptr += sizeof(int);

        memcpy(&fgm->Nvar, ptr, sizeof(int));
        ptr += sizeof(int);

        fgm->Ngrid = reinterpret_cast<int*>(const_cast<char*>(ptr));
        ptr += sizeof(int) * fgm->Ncv;

        fgm->varname = reinterpret_cast<char (*)[VAR_NAME_LENGTH]>(const_cast<char*>(ptr));
        ptr += VAR_NAME_LENGTH * fgm->Nvar;

        fgm->data = reinterpret_cast<float*>(const_cast<char*>(ptr));
        ptr += sizeof(float) * fgm->Ntotal * fgm->Nvar;

        return fgm;
    }

    void FGMTable::cleanup()
    {
        for (auto& it : fgmInstances_)
        {
            if (it.second)
            {
                freeFGM(it.second);
                it.second = nullptr;
            }
        }
        fgmInstances_.clear();

        for (auto& it : shmWins_)
        {
            if (it.second != MPI_WIN_NULL)
            {
                MPI_Win_free(&it.second);
                it.second = MPI_WIN_NULL;
            }
        }
        shmWins_.clear();

        shmBuffers_.clear();
        shmBufferSizes_.clear();
    }  
}

// ************************************************************************* //