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

Class
    Foam::combustionModels::baseFGM

Group
    grpCombustionModels

Description
    Dummy combustion model for 'no combustion'

SourceFiles
    baseFGM.C

\*---------------------------------------------------------------------------*/

#ifndef baseFGM_H
#define baseFGM_H

#include "ThermoCombustion.H"
#include "fvOptions.H"
#include "hashedWordList.H"
#include "clockTime.H"

extern "C"
{
    #include "FGMlib.h"
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{

/*---------------------------------------------------------------------------*\
                        Class baseFGM Declaration
\*---------------------------------------------------------------------------*/

template<class ReactionThermo>
class baseFGM
:
    public ThermoCombustion<ReactionThermo>
{

protected:

    // Protected data

        // Control variable related variables

            //- Mixture fraction
            volScalarField& Z_;
            
            //- Enthalpy
            volScalarField& h_;
            
            //- Progress variable
            volScalarField& PV_;
            
            //- Progress variable source term
            volScalarField& sourcePV_;

            //- Progress variable preferential diffusion contribution of H
            volScalarField& prefDPV_H_;

            //- Progress variable preferential diffusion contribution of H2
            volScalarField& prefDPV_H2_;

            //- Progress variable preferential diffusion contribution of H2O
            volScalarField& prefDPV_H2O_;

            //- Progress variable preferential diffusion contribution of thermal diffusion
            volScalarField& prefDPV_T_;

            //- Enthalpy preferential diffusion contribution of H
            volScalarField& prefDh_H_;

            //- Enthalpy preferential diffusion contribution of H2
            volScalarField& prefDh_H2_;

            //- Enthalpy preferential diffusion contribution of H2O
            volScalarField& prefDh_H2O_;

            //- Enthalpy preferential diffusion contribution of thermal diffusion
            volScalarField& prefDh_T_;

            //- Mixture fraction preferential diffusion contribution of H
            volScalarField& prefDZ_H_;

            //- Mixture fraction preferential diffusion contribution of H2
            volScalarField& prefDZ_H2_;

            //- Mixture fraction preferential diffusion contribution of H2O
            volScalarField& prefDZ_H2O_;

            //- Mixture fraction preferential diffusion contribution of thermal diffusion
            volScalarField& prefDZ_T_;

            //- Mass fraction of hydrogen radical 
            volScalarField& H_;

            //- Mass fraction of hydrogen 
            volScalarField& H2_;

            //- Mass fraction of water 
            volScalarField& H2O_;

        // Thermo variables

            //- Temperature
            volScalarField& T_;

            //- Dynamic viscosity
            volScalarField& mu_;

            //- Thermal diffusivity
            volScalarField& alpha_;

            //- Compressibility
            volScalarField& psi_;
 
        // Look-up

            //- Create pointer to the FGM look-up routine
            FGM* fgm_;

            //- Control variables
            double* controlVariables_;

            //- Total variables
            double* variables_;

    // Protected Member Functions
        
        //- Update the thermophysical parameters
        virtual void update();

public:

    //- Runtime type information
    TypeName("baseFGM");


    // Constructors

        //- Construct from components
        baseFGM
        (
            const word& modelType,
            ReactionThermo& thermo,
            const compressibleTurbulenceModel& turb,
            const word& combustionProperties
        );


    //- Destructor
    virtual ~baseFGM();


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
    #include "baseFGM.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
