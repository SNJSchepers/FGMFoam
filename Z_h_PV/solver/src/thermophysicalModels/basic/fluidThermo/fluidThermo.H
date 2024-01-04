/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2017 OpenFOAM Foundation
    Copyright (C) 2017 OpenCFD Ltd.
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
    Foam::fluidThermo

Description
    Fundamental fluid thermodynamic properties

SourceFiles
    fluidThermo.C

\*---------------------------------------------------------------------------*/

#ifndef fluidThermo_H
#define fluidThermo_H

#include "basicThermo.H"
#include "compressibleTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class fluidThermo Declaration
\*---------------------------------------------------------------------------*/

class fluidThermo
:
    public basicThermo,
    public compressibleTransportModel
{

protected:

    // - Mixture fraction [-]
    volScalarField Z_;
    
    // - Enthalpy [J/kg]
    volScalarField h_;
    
    // - Progress variable [-]
    volScalarField PV_;

    // - Progress variable source term [kg/(m^3 s)]
    volScalarField sourcePV_;

    // - Progress variable preferential diffusion contribution of H [kg/(m s)]
    volScalarField prefDPV_H_;
    
    // - Progress variable preferential diffusion contribution of H2 [kg/(m s)]
    volScalarField prefDPV_H2_;
    
    // - Progress variable preferential diffusion contribution of H2O [kg/(m s)]
    volScalarField prefDPV_H2O_;
    
    // - Progress variable preferential diffusion contribution of thermal diffusion [kg/(K m s)]
    volScalarField prefDPV_T_;
    
    // - Enthalpy preferential diffusion contribution of H [J/(m s)]
    volScalarField prefDh_H_;
    
    // - Enthalpy preferential diffusion contribution of H2 [J/(m s)]
    volScalarField prefDh_H2_;
    
    // - Enthalpy preferential diffusion contribution of H2O [J/(m s)]
    volScalarField prefDh_H2O_;
    
    // - Enthalpy preferential diffusion contribution of thermal diffusion [J/(K m s)]
    volScalarField prefDh_T_;
    
    // - Mixture fraction preferential diffusion contribution of H [kg/(m s)]
    volScalarField prefDZ_H_;
    
    // - Mixture fraction preferential diffusion contribution of H2 [kg/(m s)]
    volScalarField prefDZ_H2_;
    
    // - Mixture fraction preferential diffusion contribution of H2O [kg/(m s)]
    volScalarField prefDZ_H2O_;
    
    // - Mixture fraction preferential diffusion contribution of thermal diffusion [kg/(K m s)]
    volScalarField prefDZ_T_;
    
    // - Mass fraction of hydrogen radical [-]
    volScalarField H_;
    
    // - Mass fraction of hydrogen [-]
    volScalarField H2_;
    
    // - Mass fraction of water [-]
    volScalarField H2O_;
    
public:

    //- Runtime type information
    TypeName("fluidThermo");


    //- Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        fluidThermo,
        fvMesh,
        (const fvMesh& mesh, const word& phaseName),
        (mesh, phaseName)
    );


    //- Declare run-time constructor selection table for dictionary based
    declareRunTimeSelectionTable
    (
        autoPtr,
        fluidThermo,
        fvMeshDictPhase,
        (const fvMesh& mesh, const word& phaseName, const word& dictName),
        (mesh, phaseName, dictName)
    );

    // Constructors

        //- Construct from mesh and phase name
        fluidThermo
        (
            const fvMesh&,
            const word& phaseName
        );

        //- Construct from mesh and phase name
        fluidThermo
        (
            const fvMesh&,
            const dictionary&,
            const word& phaseName
        );


        //- Construct from mesh,dictionary,phase name with a single temperature
        fluidThermo
        (
            const fvMesh&,
            const word& phaseName,
            const word& dictionaryName
        );


    //- Selector
    static autoPtr<fluidThermo> New
    (
        const fvMesh&,
        const word& phaseName=word::null
    );


    //- Selector
    static autoPtr<fluidThermo> New
    (
        const fvMesh&,
        const word& phaseName,
        const word& dictName
    );


    //- Destructor
    virtual ~fluidThermo();


    // Member functions

        // Access to thermodynamic state variables

            //- Add the given density correction to the density field.
            //  Used to update the density field following pressure solution
            //  Limit thermo rho between rhoMin and rhoMax
            virtual void correctRho
            (
                const volScalarField& deltaRho,
                const dimensionedScalar& rhoMin,
                const dimensionedScalar& rhoMax
            ) = 0;


            //- Add the given density correction to the density field.
            //  Used to update the density field following pressure solution
            virtual void correctRho
            (
                const volScalarField& deltaRho
            ) = 0;

            //- Compressibility [s^2/m^2]
            virtual const volScalarField& psi() const = 0;


        // Access to transport state variables

            //- Dynamic viscosity of mixture [kg/m/s]
            virtual tmp<volScalarField> mu() const = 0;

            //- Dynamic viscosity of mixture for patch [kg/m/s]
            virtual tmp<scalarField> mu(const label patchi) const = 0;

            //- Kinematic viscosity of mixture [m^2/s]
            virtual tmp<volScalarField> nu() const;

            //- Kinematic viscosity of mixture for patch [m^2/s]
            virtual tmp<scalarField> nu(const label patchi) const;
       
       // Access to the control variables
        
            //- Mixture fraction [-]
            virtual volScalarField& Z();
            virtual const volScalarField& Z() const;

            //- Enthalpy [J/kg]
            virtual volScalarField& h();
            virtual const volScalarField& h() const;

            //- Progress variable [-]
            virtual volScalarField& PV();
            virtual const volScalarField& PV() const;
            
       // Access to the dependent variables

            //- Progress variable source term [kg/(m^3 s)]
            virtual volScalarField& sourcePV();
            virtual const volScalarField& sourcePV() const;

            //- Progress variable preferential diffusion contribution of H [kg/(m s)]
            virtual volScalarField& prefDPV_H();
            virtual const volScalarField& prefDPV_H() const;
            
            //- Progress variable preferential diffusion contribution of H2 [kg/(m s)]
            virtual volScalarField& prefDPV_H2();
            virtual const volScalarField& prefDPV_H2() const;
            
            //- Progress variable preferential diffusion contribution of H2O [kg/(m s)]
            virtual volScalarField& prefDPV_H2O();
            virtual const volScalarField& prefDPV_H2O() const;
            
            //- Progress variable preferential diffusion contribution of thermal diffusion [kg/(K m s)]
            virtual volScalarField& prefDPV_T();
            virtual const volScalarField& prefDPV_T() const;
            
            //- Enthalpy preferential diffusion contribution of H [J/(m s)]
            virtual volScalarField& prefDh_H();
            virtual const volScalarField& prefDh_H() const;
            
            //- Enthalpy preferential diffusion contribution of H2 [J/(m s)]
            virtual volScalarField& prefDh_H2();
            virtual const volScalarField& prefDh_H2() const;
            
            //- Enthalpy preferential diffusion contribution of H2O [J/(m s)]
            virtual volScalarField& prefDh_H2O();
            virtual const volScalarField& prefDh_H2O() const;
            
            //- Enthalpy preferential diffusion contribution of thermal diffusion [kg/(K m s)]
            virtual volScalarField& prefDh_T();
            virtual const volScalarField& prefDh_T() const;
            
            //- Mixture fraction preferential diffusion contribution of H [kg/(m s)]
            virtual volScalarField& prefDZ_H();
            virtual const volScalarField& prefDZ_H() const;
            
            //- Mixture fraction preferential diffusion contribution of H2 [kg/(m s)]
            virtual volScalarField& prefDZ_H2();
            virtual const volScalarField& prefDZ_H2() const;
            
            //- Mixture fraction preferential diffusion contribution of H2O [kg/(m s)]
            virtual volScalarField& prefDZ_H2O();
            virtual const volScalarField& prefDZ_H2O() const;
            
            //- Mixture fraction preferential diffusion contribution of thermal diffusion [kg/(K m s)]
            virtual volScalarField& prefDZ_T();
            virtual const volScalarField& prefDZ_T() const;
            
            //- Mass fraction of hydrogen radical [-]
            virtual volScalarField& H();
            virtual const volScalarField& H() const;
            
            //- Mass fraction of hydrogen [-]
            virtual volScalarField& H2();
            virtual const volScalarField& H2() const;
            
            //- Mass fraction of water [-]
            virtual volScalarField& H2O();
            virtual const volScalarField& H2O() const;
            
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
