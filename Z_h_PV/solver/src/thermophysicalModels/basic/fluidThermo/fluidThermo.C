/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "fluidThermo.H"

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(fluidThermo, 0);
    defineRunTimeSelectionTable(fluidThermo, fvMesh);
    defineRunTimeSelectionTable(fluidThermo, fvMeshDictPhase);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluidThermo::fluidThermo(const fvMesh& mesh, const word& phaseName)
:
    basicThermo(mesh, phaseName),

    Z_
    (
        IOobject
        (
            phasePropertyName("Z"),
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
    mesh
    ),

    h_
    (
        IOobject
        (
            phasePropertyName("h"),
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
    mesh
    ),

    PV_
    (
        IOobject
        (
            phasePropertyName("PV"),
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
    mesh
    ),

    sourcePV_
    (
        IOobject
        (
            phasePropertyName("sourcePV"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    mesh,
    dimensionedScalar(dimDensity/dimTime, 0)
    ),

    prefDPV_H_
    (
        IOobject
        (
            phasePropertyName("prefDPV_H"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    mesh,
    dimensionedScalar(dimMass/dimLength/dimTime, 0)
    ),

    prefDPV_H2_
    (
        IOobject
        (
            phasePropertyName("prefDPV_H2"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    mesh,
    dimensionedScalar(dimMass/dimLength/dimTime, 0)
    ),

    prefDPV_H2O_
    (
        IOobject
        (
            phasePropertyName("prefDPV_H2O"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    mesh,
    dimensionedScalar(dimMass/dimLength/dimTime, 0)
    ),

    prefDPV_T_
    (
        IOobject
        (
            phasePropertyName("prefDPV_T"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    mesh,
    dimensionedScalar(dimMass/dimLength/dimTime/dimTemperature, 0)
    ),

    prefDh_H_
    (
        IOobject
        (
            phasePropertyName("prefDh_H"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    mesh,
    dimensionedScalar(dimEnergy/dimLength/dimTime, 0)
    ),

    prefDh_H2_
    (
        IOobject
        (
            phasePropertyName("prefDh_H2"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    mesh,
    dimensionedScalar(dimEnergy/dimLength/dimTime, 0)
    ),

    prefDh_H2O_
    (
        IOobject
        (
            phasePropertyName("prefDh_H2O"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    mesh,
    dimensionedScalar(dimEnergy/dimLength/dimTime, 0)
    ),

    prefDh_T_
    (
        IOobject
        (
            phasePropertyName("prefDh_T"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    mesh,
    dimensionedScalar(dimEnergy/dimLength/dimTime/dimTemperature, 0)
    ),

    prefDZ_H_
    (
        IOobject
        (
            phasePropertyName("prefDZ_H"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    mesh,
    dimensionedScalar(dimMass/dimLength/dimTime, 0)
    ),

    prefDZ_H2_
    (
        IOobject
        (
            phasePropertyName("prefDZ_H2"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    mesh,
    dimensionedScalar(dimMass/dimLength/dimTime, 0)
    ),

    prefDZ_H2O_
    (
        IOobject
        (
            phasePropertyName("prefDZ_H2O"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    mesh,
    dimensionedScalar(dimMass/dimLength/dimTime, 0)
    ),

    prefDZ_T_
    (
        IOobject
        (
            phasePropertyName("prefDZ_T"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    mesh,
    dimensionedScalar(dimMass/dimLength/dimTime/dimTemperature, 0)
    ),

    H_
    (
        IOobject
        (
            phasePropertyName("H"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    mesh,
    dimensionedScalar(dimless, 0)
    ),

    H2_
    (
        IOobject
        (
            phasePropertyName("H2"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    mesh,
    dimensionedScalar(dimless, 0)
    ),

    H2O_
    (
        IOobject
        (
            phasePropertyName("H2O"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    mesh,
    dimensionedScalar(dimless, 0)
    )
{}



Foam::fluidThermo::fluidThermo
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName
)
:
    basicThermo(mesh, dict, phaseName),

    Z_
    (
        IOobject
        (
            phasePropertyName("Z"),
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
    mesh
    ),

    h_
    (
        IOobject
        (
            phasePropertyName("h"),
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
    mesh
    ),

    PV_
    (
        IOobject
        (
            phasePropertyName("PV"),
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
    mesh
    ),

    sourcePV_
    (
        IOobject
        (
            phasePropertyName("sourcePV"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    mesh,
    dimensionedScalar(dimDensity/dimTime, 0)
    ),

    prefDPV_H_
    (
        IOobject
        (
            phasePropertyName("prefDPV_H"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    mesh,
    dimensionedScalar(dimMass/dimLength/dimTime, 0)
    ),

    prefDPV_H2_
    (
        IOobject
        (
            phasePropertyName("prefDPV_H2"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    mesh,
    dimensionedScalar(dimMass/dimLength/dimTime, 0)
    ),

    prefDPV_H2O_
    (
        IOobject
        (
            phasePropertyName("prefDPV_H2O"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    mesh,
    dimensionedScalar(dimMass/dimLength/dimTime, 0)
    ),

    prefDPV_T_
    (
        IOobject
        (
            phasePropertyName("prefDPV_T"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    mesh,
    dimensionedScalar(dimMass/dimLength/dimTime/dimTemperature, 0)
    ),

    prefDh_H_
    (
        IOobject
        (
            phasePropertyName("prefDh_H"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    mesh,
    dimensionedScalar(dimEnergy/dimLength/dimTime, 0)
    ),

    prefDh_H2_
    (
        IOobject
        (
            phasePropertyName("prefDh_H2"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    mesh,
    dimensionedScalar(dimEnergy/dimLength/dimTime, 0)
    ),

    prefDh_H2O_
    (
        IOobject
        (
            phasePropertyName("prefDh_H2O"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    mesh,
    dimensionedScalar(dimEnergy/dimLength/dimTime, 0)
    ),

    prefDh_T_
    (
        IOobject
        (
            phasePropertyName("prefDh_T"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    mesh,
    dimensionedScalar(dimEnergy/dimLength/dimTime/dimTemperature, 0)
    ),

    prefDZ_H_
    (
        IOobject
        (
            phasePropertyName("prefDZ_H"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    mesh,
    dimensionedScalar(dimMass/dimLength/dimTime, 0)
    ),

    prefDZ_H2_
    (
        IOobject
        (
            phasePropertyName("prefDZ_H2"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    mesh,
    dimensionedScalar(dimMass/dimLength/dimTime, 0)
    ),

    prefDZ_H2O_
    (
        IOobject
        (
            phasePropertyName("prefDZ_H2O"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    mesh,
    dimensionedScalar(dimMass/dimLength/dimTime, 0)
    ),

    prefDZ_T_
    (
        IOobject
        (
            phasePropertyName("prefDZ_T"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    mesh,
    dimensionedScalar(dimMass/dimLength/dimTime/dimTemperature, 0)
    ),

    H_
    (
        IOobject
        (
            phasePropertyName("H"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    mesh,
    dimensionedScalar(dimless, 0)
    ),

    H2_
    (
        IOobject
        (
            phasePropertyName("H2"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    mesh,
    dimensionedScalar(dimless, 0)
    ),

    H2O_
    (
        IOobject
        (
            phasePropertyName("H2O"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    mesh,
    dimensionedScalar(dimless, 0)
    )
{}


Foam::fluidThermo::fluidThermo
(
    const fvMesh& mesh,
    const word& phaseName,
    const word& dictionaryName
)
:
    basicThermo(mesh, phaseName, dictionaryName),

    Z_
    (
        IOobject
        (
            phasePropertyName("Z"),
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
    mesh
    ),

    h_
    (
        IOobject
        (
            phasePropertyName("h"),
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
    mesh
    ),

    PV_
    (
        IOobject
        (
            phasePropertyName("PV"),
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
    mesh
    ),

    sourcePV_
    (
        IOobject
        (
            phasePropertyName("sourcePV"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    mesh,
    dimensionedScalar(dimDensity/dimTime, 0)
    ),

    prefDPV_H_
    (
        IOobject
        (
            phasePropertyName("prefDPV_H"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    mesh,
    dimensionedScalar(dimMass/dimLength/dimTime, 0)
    ),

    prefDPV_H2_
    (
        IOobject
        (
            phasePropertyName("prefDPV_H2"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    mesh,
    dimensionedScalar(dimMass/dimLength/dimTime, 0)
    ),

    prefDPV_H2O_
    (
        IOobject
        (
            phasePropertyName("prefDPV_H2O"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    mesh,
    dimensionedScalar(dimMass/dimLength/dimTime, 0)
    ),

    prefDPV_T_
    (
        IOobject
        (
            phasePropertyName("prefDPV_T"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    mesh,
    dimensionedScalar(dimMass/dimLength/dimTime/dimTemperature, 0)
    ),

    prefDh_H_
    (
        IOobject
        (
            phasePropertyName("prefDh_H"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    mesh,
    dimensionedScalar(dimEnergy/dimLength/dimTime, 0)
    ),

    prefDh_H2_
    (
        IOobject
        (
            phasePropertyName("prefDh_H2"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    mesh,
    dimensionedScalar(dimEnergy/dimLength/dimTime, 0)
    ),

    prefDh_H2O_
    (
        IOobject
        (
            phasePropertyName("prefDh_H2O"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    mesh,
    dimensionedScalar(dimEnergy/dimLength/dimTime, 0)
    ),

    prefDh_T_
    (
        IOobject
        (
            phasePropertyName("prefDh_T"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    mesh,
    dimensionedScalar(dimEnergy/dimLength/dimTime/dimTemperature, 0)
    ),

    prefDZ_H_
    (
        IOobject
        (
            phasePropertyName("prefDZ_H"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    mesh,
    dimensionedScalar(dimMass/dimLength/dimTime, 0)
    ),

    prefDZ_H2_
    (
        IOobject
        (
            phasePropertyName("prefDZ_H2"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    mesh,
    dimensionedScalar(dimMass/dimLength/dimTime, 0)
    ),

    prefDZ_H2O_
    (
        IOobject
        (
            phasePropertyName("prefDZ_H2O"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    mesh,
    dimensionedScalar(dimMass/dimLength/dimTime, 0)
    ),

    prefDZ_T_
    (
        IOobject
        (
            phasePropertyName("prefDZ_T"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    mesh,
    dimensionedScalar(dimMass/dimLength/dimTime/dimTemperature, 0)
    ),

    H_
    (
        IOobject
        (
            phasePropertyName("H"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    mesh,
    dimensionedScalar(dimless, 0)
    ),

    H2_
    (
        IOobject
        (
            phasePropertyName("H2"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    mesh,
    dimensionedScalar(dimless, 0)
    ),

    H2O_
    (
        IOobject
        (
            phasePropertyName("H2O"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    mesh,
    dimensionedScalar(dimless, 0)
    )
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::fluidThermo> Foam::fluidThermo::New
(
    const fvMesh& mesh,
    const word& phaseName
)
{
    return basicThermo::New<fluidThermo>(mesh, phaseName);
}


Foam::autoPtr<Foam::fluidThermo> Foam::fluidThermo::New
(
    const fvMesh& mesh,
    const word& phaseName,
    const word& dictName
)
{
    return basicThermo::New<fluidThermo>(mesh, phaseName, dictName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluidThermo::~fluidThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::fluidThermo::nu() const
{
    return mu()/rho();
}


Foam::tmp<Foam::scalarField> Foam::fluidThermo::nu(const label patchi) const
{
    return mu(patchi)/rho(patchi);
}

Foam::volScalarField& Foam::fluidThermo::Z()
{
    return Z_;
}

const Foam::volScalarField& Foam::fluidThermo::Z() const
{
    return Z_;
}

Foam::volScalarField& Foam::fluidThermo::h()
{
    return h_;
}

const Foam::volScalarField& Foam::fluidThermo::h() const
{
    return h_;
}

Foam::volScalarField& Foam::fluidThermo::PV()
{
    return PV_;
}

const Foam::volScalarField& Foam::fluidThermo::PV() const
{
    return PV_;
}

Foam::volScalarField& Foam::fluidThermo::sourcePV()
{
    return sourcePV_;
}

const Foam::volScalarField& Foam::fluidThermo::sourcePV() const
{
    return sourcePV_;
}

Foam::volScalarField& Foam::fluidThermo::prefDPV_H()
{
    return prefDPV_H_;
}

const Foam::volScalarField& Foam::fluidThermo::prefDPV_H() const
{
    return prefDPV_H_;
}

Foam::volScalarField& Foam::fluidThermo::prefDPV_H2()
{
    return prefDPV_H2_;
}

const Foam::volScalarField& Foam::fluidThermo::prefDPV_H2() const
{
    return prefDPV_H2_;
}

Foam::volScalarField& Foam::fluidThermo::prefDPV_H2O()
{
    return prefDPV_H2O_;
}

const Foam::volScalarField& Foam::fluidThermo::prefDPV_H2O() const
{
    return prefDPV_H2O_;
}

Foam::volScalarField& Foam::fluidThermo::prefDPV_T()
{
    return prefDPV_T_;
}

const Foam::volScalarField& Foam::fluidThermo::prefDPV_T() const
{
    return prefDPV_T_;
}

Foam::volScalarField& Foam::fluidThermo::prefDh_H()
{
    return prefDh_H_;
}

const Foam::volScalarField& Foam::fluidThermo::prefDh_H() const
{
    return prefDh_H_;
}

Foam::volScalarField& Foam::fluidThermo::prefDh_H2()
{
    return prefDh_H2_;
}

const Foam::volScalarField& Foam::fluidThermo::prefDh_H2() const
{
    return prefDh_H2_;
}

Foam::volScalarField& Foam::fluidThermo::prefDh_H2O()
{
    return prefDh_H2O_;
}

const Foam::volScalarField& Foam::fluidThermo::prefDh_H2O() const
{
    return prefDh_H2O_;
}

Foam::volScalarField& Foam::fluidThermo::prefDh_T()
{
    return prefDh_T_;
}

const Foam::volScalarField& Foam::fluidThermo::prefDh_T() const
{
    return prefDh_T_;
}

Foam::volScalarField& Foam::fluidThermo::prefDZ_H()
{
    return prefDZ_H_;
}

const Foam::volScalarField& Foam::fluidThermo::prefDZ_H() const
{
    return prefDZ_H_;
}

Foam::volScalarField& Foam::fluidThermo::prefDZ_H2()
{
    return prefDZ_H2_;
}

const Foam::volScalarField& Foam::fluidThermo::prefDZ_H2() const
{
    return prefDZ_H2_;
}

Foam::volScalarField& Foam::fluidThermo::prefDZ_H2O()
{
    return prefDZ_H2O_;
}

const Foam::volScalarField& Foam::fluidThermo::prefDZ_H2O() const
{
    return prefDZ_H2O_;
}

Foam::volScalarField& Foam::fluidThermo::prefDZ_T()
{
    return prefDZ_T_;
}

const Foam::volScalarField& Foam::fluidThermo::prefDZ_T() const
{
    return prefDZ_T_;
}

Foam::volScalarField& Foam::fluidThermo::H()
{
    return H_;
}

const Foam::volScalarField& Foam::fluidThermo::H() const
{
    return H_;
}

Foam::volScalarField& Foam::fluidThermo::H2()
{
    return H2_;
}

const Foam::volScalarField& Foam::fluidThermo::H2() const
{
    return H2_;
}

Foam::volScalarField& Foam::fluidThermo::H2O()
{
    return H2O_;
}

const Foam::volScalarField& Foam::fluidThermo::H2O() const
{
    return H2O_;
}

// ************************************************************************* //
