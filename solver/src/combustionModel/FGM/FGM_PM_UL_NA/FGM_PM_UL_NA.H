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

Class
    Foam::combustionModels::FGM_PM_UL_NA

Group
    grpCombustionModels

Description
    Dummy combustion model for 'no combustion'

SourceFiles
    FGM_PM_UL_NA.C

\*---------------------------------------------------------------------------*/

#ifndef FGM_PM_UL_NA_H
#define FGM_PM_UL_NA_H

#include "ThermoCombustion.H"
#include "fvOptions.H"
#include "hashedWordList.H"
#include "clockTime.H"
#include "wallFvPatch.H"

#include "FGMTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{

/*---------------------------------------------------------------------------*\
                        Class FGM_PM_UL_NA Declaration
\*---------------------------------------------------------------------------*/

template<class ReactionThermo>
class FGM_PM_UL_NA
:
    public ThermoCombustion<ReactionThermo>
{

protected:

    // Protected data

        // FGM fields
            
			//- Progress variable
            volScalarField Yc_;
            
            //- Enthalpy
            volScalarField ht_;
            
			//- Progress variable source term
            volScalarField sourceYc_;
			
			//- Methane mass fraction
            volScalarField CH4_;
			
            //- Hydrogen mass fraction
            volScalarField H2_;
            
            //- Water mass fraction
            volScalarField H2O_;

            //- Hydroxyl mass fraction
            volScalarField OH_;

            //- Oxygen mass fraction
            volScalarField O2_;

            //- Carbon dioxide mass fraction
            volScalarField CO2_;

            //- Carbon monoxide mass fraction
            volScalarField CO_;

            //- Heat release rate
            volScalarField HeatRelease_;
 
        // Look-up

            //- Create pointer to the FGM look-up routine
            FGM* fgm_;
            
            //- Control variables
            float* controlVariables_;

            //- Total variables
            float* variables_;

            //- Indices for the dependent variables in the FGM table
            int isourceYc; 
            int iCH4;
            int iH2;
            int iH2O;
            int iOH;
            int iO2;
            int iCO2;
            int iCO;
            int iT;
            int imu;
            int icp;
            int ilambda;
            int iHeatRelease;
            int irho;

        // Other

            //- Turbulent Schmidt Number
            scalar Sct_;

    // Protected Member Functions

        //- solve the control variable transport equations
        virtual void solve();
    
        //- Update the thermophysical parameters
        virtual void update();


public:

    //- Runtime type information
    TypeName("FGM_PM_UL_NA");


    // Constructors

        //- Construct from components
        FGM_PM_UL_NA
        (
            const word& modelType,
            ReactionThermo& thermo,
            const compressibleTurbulenceModel& turb,
            const word& combustionProperties
        );


    //- Destructor
    virtual ~FGM_PM_UL_NA();


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
    #include "FGM_PM_UL_NA.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
