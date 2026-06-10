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
    Copyright (C) 2019 OpenCFD Ltd.
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

#include "FGM_PM_UL_NA.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ReactionThermo>
Foam::combustionModels::FGM_PM_UL_NA<ReactionThermo>::FGM_PM_UL_NA
(
    const word& modelType,
    ReactionThermo& thermo,
    const compressibleTurbulenceModel& turb,
    const word& combustionProperties
)
:
    ThermoCombustion<ReactionThermo>(modelType, thermo, turb),
    Yc_
    (
        IOobject
        (
            "Yc",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
    this->mesh()
    ),
    ht_
    (
        IOobject
        (
            "ht",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh()
    ),
    sourceYc_
    (
        IOobject
        (
            "sourceYc",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        dimensionedScalar("zero", dimDensity/dimTime, 0.0)
    ),
    CH4_
    (
        IOobject
        (
            "CH4",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        dimensionedScalar("zero", dimless, 0.0)
    ),
    H2_
    (
        IOobject
        (
            "H2",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        dimensionedScalar("zero", dimless, 0.0)
    ),
    H2O_
    (
        IOobject
        (
            "H2O",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        dimensionedScalar("zero", dimless, 0.0)
    ),
    OH_
    (
        IOobject
        (
            "OH",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        dimensionedScalar("zero", dimless, 0.0)
    ),
    O2_
    (
        IOobject
        (
            "O2",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        dimensionedScalar("zero", dimless, 0.0)
    ),
    CO2_
    (
        IOobject
        (
            "CO2",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        dimensionedScalar("zero", dimless, 0.0)
    ),
    CO_
    (
        IOobject
        (
            "CO",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        dimensionedScalar("zero", dimless, 0.0)
    ),
    HeatRelease_
    (
        IOobject
        (
            "HeatRelease",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        dimensionedScalar("zero", dimEnergy/dimVolume/dimTime, 0.0)
    ),

	fgm_(nullptr),

    isourceYc(-1),
    iCH4(-1),
    iH2(-1),
    iH2O(-1),
    iOH(-1),
    iO2(-1),
    iCO2(-1),
    iCO(-1),
    iT(-1),
    imu(-1),
    icp(-1),
    ilambda(-1),
    iHeatRelease(-1),
    irho(-1)
{
    // Initialise control and total variables as null pointer (empty)
    controlVariables_ = nullptr;
    variables_        = nullptr;
		
    // Release memory if arrays already exist
    if (controlVariables_) {
        delete[] controlVariables_;
    }
    if (variables_) {
        delete[] variables_;
    }

    fgm_ = FGMTable::getFGMInstance("constant/lookUp/database.fgm");
	
	// Initialize arrays with sizes obtained from the FGM table
    controlVariables_ = new float[fgm_->Ncv];
    variables_        = new float[fgm_->Nvar];

    for (int i = fgm_->Ncv; i < fgm_->Nvar; i++)
    {
        if (std::strcmp(fgm_->varname[i], "SOURCE_CV1") == 0){
            isourceYc = i; 
        }
        if (std::strcmp(fgm_->varname[i], "CH4") == 0){
            iCH4 = i;
        }
        if (std::strcmp(fgm_->varname[i], "H2") == 0){
            iH2 = i;
        }
        if (std::strcmp(fgm_->varname[i], "H2O") == 0){
            iH2O = i;
        }
        if (std::strcmp(fgm_->varname[i], "OH") == 0){
            iOH = i;
        }
        if (std::strcmp(fgm_->varname[i], "O2") == 0){
            iO2 = i;
        }
        if (std::strcmp(fgm_->varname[i], "CO2") == 0){
            iCO2 = i;
        }
        if (std::strcmp(fgm_->varname[i], "CO") == 0){
            iCO = i;
        }
        if (std::strcmp(fgm_->varname[i], "TEMPERATURE") == 0){
            iT = i;
        }
        if (std::strcmp(fgm_->varname[i], "VISCOSITY") == 0){
            imu = i;
        }
        if (std::strcmp(fgm_->varname[i], "CP") == 0){
            icp = i;
        }
        if (std::strcmp(fgm_->varname[i], "CONDUCTIVITY") == 0){
            ilambda = i;
        }
        if (std::strcmp(fgm_->varname[i], "HEATRELEASE") == 0){
            iHeatRelease = i;
        }
        if (std::strcmp(fgm_->varname[i], "DENSITY") == 0){
            irho = i;
        }
    }
	
    // Update the field variables with the table variables
    update();

    // Read turbulent Schmidt number from combustionProperties if present, otherwise default to 0.7
    if (this->coeffs().found("turbulenceCoeffs"))
    {
        Sct_ = this->coeffs().subDict("turbulenceCoeffs")
                   .template lookupOrDefault<scalar>("Sct", 0.7);
    }
    else
    {
        Sct_ = 0.7;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ReactionThermo>
Foam::combustionModels::FGM_PM_UL_NA<ReactionThermo>::~FGM_PM_UL_NA()
{
    FGMTable::cleanup();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ReactionThermo>
void Foam::combustionModels::FGM_PM_UL_NA<ReactionThermo>::solve()
{
    
    // Access to fields
    const volScalarField& rho     = this->mesh().template lookupObject<volScalarField>("rho");
    const volScalarField& alpha   = this->mesh().template lookupObject<volScalarField>("alpha");
    const surfaceScalarField& phi = this->mesh().template lookupObject<surfaceScalarField>("phi");

    // Field table for multivariate interpolation of control variables
    multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

    // Add fields to the list for convection scheme
    fields.add(Yc_);
    fields.add(ht_);

    tmp<fv::convectionScheme<scalar>> mvConvection =
        fv::convectionScheme<scalar>::New
        (
            this->mesh(),
            fields,
            phi,
            this->mesh().divScheme("div(phi,Yc)")
        );

    // Solve the progress variable equation
    fvScalarMatrix YcEqn
    (
        fvm::ddt(rho, Yc_)
      + mvConvection().fvmDiv(phi, Yc_)
      - fvm::laplacian(alpha + this->turbulence().mut()/Sct_, Yc_)
     ==
        sourceYc_
    );
    YcEqn.relax();
    YcEqn.solve();

    // Solve the enthalpy equation
    fvScalarMatrix htEqn
    (
        fvm::ddt(rho, ht_)
      + mvConvection().fvmDiv(phi, ht_)
      - fvm::laplacian(alpha + this->turbulence().mut()/Sct_, ht_)
    );
    htEqn.relax();
    htEqn.solve();
}

template<class ReactionThermo>
void Foam::combustionModels::FGM_PM_UL_NA<ReactionThermo>::update()
{
    const clockTime clockTime_= clockTime();
    clockTime_.timeIncrement();

    // Create reference to the thermo variables
    volScalarField& T     = const_cast<volScalarField&> ( this->mesh().template lookupObject<volScalarField>("T") );
    volScalarField& psi   = const_cast<volScalarField&> ( this->mesh().template lookupObject<volScalarField>("thermo:psi") );
    volScalarField& mu    = const_cast<volScalarField&> ( this->mesh().template lookupObject<volScalarField>("mu") );
    volScalarField& alpha = const_cast<volScalarField&> ( this->mesh().template lookupObject<volScalarField>("alpha") );

    // Create reference to the field cell center values
    scalarField& YcCells          = this->Yc_.primitiveFieldRef();
    scalarField& htCells          = this->ht_.primitiveFieldRef();

	scalarField& sourceYcCells    = this->sourceYc_.primitiveFieldRef();
	scalarField& CH4Cells         = this->CH4_.primitiveFieldRef();
	scalarField& H2Cells          = this->H2_.primitiveFieldRef();
	scalarField& H2OCells         = this->H2O_.primitiveFieldRef();
	scalarField& OHCells          = this->OH_.primitiveFieldRef();
	scalarField& O2Cells          = this->O2_.primitiveFieldRef();
	scalarField& CO2Cells         = this->CO2_.primitiveFieldRef();
	scalarField& COCells          = this->CO_.primitiveFieldRef();
	scalarField& HeatReleaseCells = this->HeatRelease_.primitiveFieldRef();
    scalarField& TCells           = T.primitiveFieldRef();
    scalarField& muCells          = mu.primitiveFieldRef();
    scalarField& alphaCells       = alpha.primitiveFieldRef();
    scalarField& psiCells         = psi.primitiveFieldRef();

    // Loop over all cells
    forAll(YcCells, celli)
    {
        // Scale the control variables
        controlVariables_[0] = YcCells[celli];
        controlVariables_[1] = htCells[celli];

        if (lookupFGM_2D(fgm_,controlVariables_,variables_) == EXIT_SUCCESS)
        {
            sourceYcCells[celli]    = variables_[isourceYc]; 
            CH4Cells[celli]         = variables_[iCH4];
            H2Cells[celli]          = variables_[iH2];
            H2OCells[celli]         = variables_[iH2O];
            OHCells[celli]          = variables_[iOH];
            O2Cells[celli]          = variables_[iO2];
            CO2Cells[celli]         = variables_[iCO2];
            COCells[celli]          = variables_[iCO];
            TCells[celli]           = variables_[iT];
            muCells[celli]          = variables_[imu];
            alphaCells[celli]       = variables_[ilambda]/variables_[icp];
            HeatReleaseCells[celli] = variables_[iHeatRelease];
            psiCells[celli]         = variables_[irho]/fgm_->pressure;
            
            // Source term of the progress variable is strictly positive
            if (sourceYcCells[celli] < 0.0)
            {
                sourceYcCells[celli]   = 0.0;
            }
        }
    }

    forAll(Yc_.boundaryFieldRef(), patchi)
    {
        fvPatchScalarField& pYc          = this->Yc_.boundaryFieldRef()[patchi];
        fvPatchScalarField& pht          = this->ht_.boundaryFieldRef()[patchi];

        fvPatchScalarField& psourceYc    = this->sourceYc_.boundaryFieldRef()[patchi];
        fvPatchScalarField& pCH4         = this->CH4_.boundaryFieldRef()[patchi];
        fvPatchScalarField& pH2          = this->H2_.boundaryFieldRef()[patchi];
        fvPatchScalarField& pH2O         = this->H2O_.boundaryFieldRef()[patchi];
        fvPatchScalarField& pOH          = this->OH_.boundaryFieldRef()[patchi];
        fvPatchScalarField& pO2          = this->O2_.boundaryFieldRef()[patchi];
        fvPatchScalarField& pCO2         = this->CO2_.boundaryFieldRef()[patchi];
        fvPatchScalarField& pCO          = this->CO_.boundaryFieldRef()[patchi];
        fvPatchScalarField& pHeatRelease = this->HeatRelease_.boundaryFieldRef()[patchi];
        fvPatchScalarField& pT           = T.boundaryFieldRef()[patchi];
        fvPatchScalarField& pmu          = mu.boundaryFieldRef()[patchi];
        fvPatchScalarField& palpha       = alpha.boundaryFieldRef()[patchi];
        fvPatchScalarField& ppsi         = psi.boundaryFieldRef()[patchi];

        forAll(pT, facei)
        {
            controlVariables_[0] = pYc[facei];
            controlVariables_[1] = pht[facei];
            
            if (lookupFGM_2D(fgm_,controlVariables_,variables_) == EXIT_SUCCESS)
            {
                psourceYc[facei]    = variables_[isourceYc];
                pCH4[facei]         = variables_[iCH4];
                pH2[facei]          = variables_[iH2];
                pH2O[facei]         = variables_[iH2O];
                pOH[facei]          = variables_[iOH];
                pO2[facei]          = variables_[iO2];
                pCO2[facei]         = variables_[iCO2];
                pCO[facei]          = variables_[iCO];
                pmu[facei]          = variables_[imu];
                palpha[facei]       = variables_[ilambda]/variables_[icp];
                pHeatRelease[facei] = variables_[iHeatRelease];
                ppsi[facei]         = variables_[irho]/fgm_->pressure;

                if (!pT.fixesValue()){
                    pT[facei] = variables_[iT];
                }

                // No reaction at the boundary faces
                if(this->mesh().boundary()[patchi].type() == "wall"){
                    psourceYc[facei]  = 0.0;
                }

                // Source term of the progress variable is strictly positive
                if (psourceYc[facei] < 0.0){
                    psourceYc[facei]  = 0.0;
                }
            }
        }
    }
    
    Info << "Parameter update time: (" << clockTime_.timeIncrement() << " s)" << endl;
}


template<class ReactionThermo>
void Foam::combustionModels::FGM_PM_UL_NA<ReactionThermo>::correct()
{
    solve();
    update();
}

template<class ReactionThermo>
Foam::tmp<Foam::fvScalarMatrix>
Foam::combustionModels::FGM_PM_UL_NA<ReactionThermo>::R
(
    volScalarField& Y
) const
{
    tmp<fvScalarMatrix> tSu
    (
        new fvScalarMatrix(Y, dimMass/dimTime)
    );

    return tSu;
}


template<class ReactionThermo>
Foam::tmp<Foam::volScalarField>
Foam::combustionModels::FGM_PM_UL_NA<ReactionThermo>::Qdot() const
{
    return tmp<volScalarField>::New
    (
        IOobject
        (
            this->thermo().phasePropertyName(typeName + ":Qdot"),
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        this->mesh(),
        dimensionedScalar(dimEnergy/dimVolume/dimTime, Zero)
    );
}

// ************************************************************************* //
