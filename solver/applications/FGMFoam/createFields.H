#include "createRDeltaT.H"


// Thermo model
Info<< "Reading thermophysical properties\n" << endl;
autoPtr<psiReactionThermo> pThermo(psiReactionThermo::New(mesh));
psiReactionThermo& thermo = pThermo();
thermo.validate(args.executable(), "h", "e");

// Pressure
volScalarField& p = thermo.p();

// Density
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    thermo.rho()
);

// Viscosity
volScalarField mu
(
    IOobject
    (
        "mu",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    thermo.mu()
);

// Thermal diffusivity
volScalarField alpha
(
    IOobject
    (
        "alpha",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    thermo.alpha()
);

// Velocity
Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


#include "compressibleCreatePhi.H"

pressureControl pressureControl(p, rho, pimple.dict(), false);
mesh.setFluxRequired(p.name());

// Turbulence model
Info<< "Creating turbulence model\n" << endl;
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);

if (turbulence->type() != "Stokes"){ // If not laminar
    Info << "WARNING: No Sub-grid scale model incorporated in this solver" << endl;
}

#include "createDpdt.H"

#include "createK.H"

#include "createMRF.H"

const dimensionedScalar rhoMax("rhoMax", dimDensity, GREAT, pimple.dict());
const dimensionedScalar rhoMin("rhoMin", dimDensity, Zero, pimple.dict());

#include "createFvOptions.H"

// Combustion model
Info<< "Creating combustion model\n" << endl;
autoPtr<CombustionModel<psiReactionThermo>> combustion
(
    CombustionModel<psiReactionThermo>::New(thermo, turbulence())
);

const word combModelName(combustion->type());
Info<< "Combustion Model Used: " << combModelName << endl;

// const volScalarField& Yc    = mesh.lookupObject<volScalarField>("Yc");

#include <fstream>
#include <string>
#include <sstream>

// std::ifstream fgmStream("constant/lookUp/database.fgm");
// if (!fgmStream.is_open())
// {
//     FatalErrorInFunction << "Could not open constant/lookUp/database.fgm" << exit(FatalError);
// }

//#include "compressibleCreatePhi.H"

thermo.T().correctBoundaryConditions();

// Update density with table value
rho = thermo.rho();

// Update mass flux with updated density
phi = linearInterpolate(rho*U) & mesh.Sf();
