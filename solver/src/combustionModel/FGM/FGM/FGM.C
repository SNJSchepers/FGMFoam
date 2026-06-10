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

\*---------------------------------------------------------------------------*/

#include "FGM.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ReactionThermo>
Foam::combustionModels::FGM<ReactionThermo>::FGM
(
    const word& modelType,
    ReactionThermo& thermo,
    const compressibleTurbulenceModel& turb,
    const word& combustionProperties
)
:
    ThermoCombustion<ReactionThermo>(modelType, thermo, turb),
    fgmmodel_(nullptr)
{
    const dictionary& coeffs = this->coeffs();

    const word flameType(coeffs.lookup("flameType"));
    const word transportModel(coeffs.lookup("transportModel"));
    const word adiabatic(coeffs.lookup("adiabatic"));
    const word subGridModel(coeffs.lookup("subGridModel"));

    if (flameType == "premixed") {
        if (transportModel == "unityLewis") {
            if (adiabatic == "no") {
                if (subGridModel == "none") {
                    fgmmodel_.reset
                    (
                        new FGM_PM_UL_NA<ReactionThermo>
                        (
                            modelType,
                            thermo,
                            turb,
                            combustionProperties
                        )
                    );
                } else {
                    FatalErrorInFunction
                        << "Unknown subGridModel option: " << subGridModel << nl
                        << "Valid options are: none"
                        << exit(FatalError);
                }
            } else {
                FatalErrorInFunction
                    << "Unknown adiabatic option: " << adiabatic << nl
                    << "Valid options are: no"
                    << exit(FatalError);
            }
        } else {
            FatalErrorInFunction
                << "Unknown transportModel option: " << transportModel << nl
                << "Valid options are: unityLewis"
                << exit(FatalError);
        }
    } else {
        FatalErrorInFunction
            << "Unknown flameType: " << flameType << nl
            << "Valid options are: premixed"
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ReactionThermo>
Foam::combustionModels::FGM<ReactionThermo>::~FGM()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * * * //

template<class ReactionThermo>
void Foam::combustionModels::FGM<ReactionThermo>::correct()
{
    fgmmodel_->correct();
}

template<class ReactionThermo>
Foam::tmp<Foam::fvScalarMatrix>
Foam::combustionModels::FGM<ReactionThermo>::R
(
    volScalarField& Y
) const
{
    return fgmmodel_->R(Y);
}


template<class ReactionThermo>
Foam::tmp<Foam::volScalarField>
Foam::combustionModels::FGM<ReactionThermo>::Qdot() const
{
    return fgmmodel_->Qdot();
}

// ************************************************************************* //

