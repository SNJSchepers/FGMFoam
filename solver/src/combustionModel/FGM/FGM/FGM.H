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
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2018 OpenFOAM Foundation
     \\/     M anipulation  |
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

Class
    Foam::FGM

Description
    FGM combustion model

SourceFiles
    FGM.C

\*---------------------------------------------------------------------------*/

#ifndef FGM_H
#define FGM_H

#include "FGM_PM_CH4_HL.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{

/*---------------------------------------------------------------------------*\
                    class FGM Declaration
\*---------------------------------------------------------------------------*/

template<class ReactionThermo>
class FGM
:
    public ThermoCombustion<ReactionThermo>
{

    // Private Data

    //- Selected internal FGM model
    autoPtr<ThermoCombustion<ReactionThermo>> fgmmodel_;

public:

    //- Runtime type information
    TypeName("FGM");

    // Constructors

        //- Construct from components
        FGM
        (
            const word& modelType,
            ReactionThermo& thermo,
            const compressibleTurbulenceModel& turb,
            const word& combustionProperties
        );

    //- Destructor
    virtual ~FGM();

    // Member Functions

        //- Correct the variables by solving equations and updating thermophysical parameters
        virtual void correct();
             
        //- Fuel consumption rate matrix
        virtual tmp<fvScalarMatrix> R(volScalarField& Y) const;

        //- Heat release rate [kg/m/s3]
        virtual tmp<volScalarField> Qdot() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace combustionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "FGM.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
