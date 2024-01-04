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

#include "baseFGM.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ReactionThermo>
Foam::combustionModels::baseFGM<ReactionThermo>::baseFGM
(
    const word& modelType,
    ReactionThermo& thermo,
    const compressibleTurbulenceModel& turb,
    const word& combustionProperties
)
:
    ThermoCombustion<ReactionThermo>(modelType, thermo, turb),
    Z_           (this->thermo().Z()),
    h_           (this->thermo().h()),
    PV_          (this->thermo().PV()),
    sourcePV_    (this->thermo().sourcePV()),
    prefDPV_H_   (this->thermo().prefDPV_H()),
    prefDPV_H2_  (this->thermo().prefDPV_H2()),
    prefDPV_H2O_ (this->thermo().prefDPV_H2O()),
    prefDPV_T_   (this->thermo().prefDPV_T()),
    prefDh_H_    (this->thermo().prefDh_H()),
    prefDh_H2_   (this->thermo().prefDh_H2()),
    prefDh_H2O_  (this->thermo().prefDh_H2O()),
    prefDh_T_    (this->thermo().prefDh_T()),
    prefDZ_H_    (this->thermo().prefDZ_H()),
    prefDZ_H2_   (this->thermo().prefDZ_H2()),
    prefDZ_H2O_  (this->thermo().prefDZ_H2O()),
    prefDZ_T_    (this->thermo().prefDZ_T()),
    H_           (this->thermo().H()),
    H2_          (this->thermo().H2()),
    H2O_         (this->thermo().H2O()),
    T_           (this->thermo().T()),
    mu_          (const_cast<volScalarField&>(this->thermo().mu()())),
    alpha_       (const_cast<volScalarField&>(this->thermo().alpha())),
    psi_         (const_cast<volScalarField&>(this->thermo().psi())),
      
    fgm_(readFGM("constant/lookUp/database.fgm"))
{
    // Initialise control and total variables as null pointer (empty)
    controlVariables_ = nullptr;
    variables_        = nullptr;

    // Release memory if arrays already exist
    if (controlVariables_) {
        delete[] controlVariables_;
    }
    if (controlVariables_) {
        delete[] controlVariables_;
    }

    // Initialize arrays with sizes obtained from the FGM table
    controlVariables_ = new double[fgm_->Ncv];
    variables_        = new double[fgm_->Nvar];

    // Update the field variables with the table variables
    update();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ReactionThermo>
Foam::combustionModels::baseFGM<ReactionThermo>::~baseFGM()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ReactionThermo>
void Foam::combustionModels::baseFGM<ReactionThermo>::update()
{
    const clockTime clockTime_= clockTime();
    clockTime_.timeIncrement();
    
    // Initialise cp and lambda at arbitrary small values. These will be overwritten by the table data 
    scalar cp     = 1.0e-10;
    scalar lambda = 1.0e-10;

    // Create reference to the field cell center values
    scalarField& ZCells           = this->Z_.primitiveFieldRef();
    scalarField& hCells           = this->h_.primitiveFieldRef();
    scalarField& PVCells          = this->PV_.primitiveFieldRef();

    scalarField& sourcePVCells    = this->sourcePV_.primitiveFieldRef();
    scalarField& prefDPV_HCells   = this->prefDPV_H_.primitiveFieldRef();
    scalarField& prefDPV_H2Cells  = this->prefDPV_H2_.primitiveFieldRef();
    scalarField& prefDPV_H2OCells = this->prefDPV_H2O_.primitiveFieldRef();
    scalarField& prefDPV_TCells   = this->prefDPV_T_.primitiveFieldRef();
    scalarField& prefDh_HCells    = this->prefDh_H_.primitiveFieldRef();
    scalarField& prefDh_H2Cells   = this->prefDh_H2_.primitiveFieldRef();
    scalarField& prefDh_H2OCells  = this->prefDh_H2O_.primitiveFieldRef();
    scalarField& prefDh_TCells    = this->prefDh_T_.primitiveFieldRef();
    scalarField& prefDZ_HCells    = this->prefDZ_H_.primitiveFieldRef();
    scalarField& prefDZ_H2Cells   = this->prefDZ_H2_.primitiveFieldRef();
    scalarField& prefDZ_H2OCells  = this->prefDZ_H2O_.primitiveFieldRef();
    scalarField& prefDZ_TCells    = this->prefDZ_T_.primitiveFieldRef();
    scalarField& HCells           = this->H_.primitiveFieldRef();
    scalarField& H2Cells          = this->H2_.primitiveFieldRef();
    scalarField& H2OCells         = this->H2O_.primitiveFieldRef();
    scalarField& TCells           = this->T_.primitiveFieldRef();
    scalarField& muCells          = this->mu_.primitiveFieldRef();
    scalarField& alphaCells       = this->alpha_.primitiveFieldRef();
    scalarField& psiCells         = this->psi_.primitiveFieldRef();

    // Loop over all cells
    forAll(PVCells, celli)
    {
        // Scale the control variables
        //controlVariables_[0] = ZCells[celli];
        //controlVariables_[1] = hCells[celli];
        //controlVariables_[2] = PVCells[celli];
        controlVariables_[0] = PVCells[celli];
        controlVariables_[1] = hCells[celli];
        controlVariables_[2] = ZCells[celli];
        
        //Info << "ZCells :"  << ZCells[celli]  << endl;
        //Info << "hCells :"  << hCells[celli]  << endl;
        //Info << "PVCells :" << PVCells[celli] << endl;

        //lookupFGM_ND(fgm_,controlVariables_,variables_);
        lookupFGM_3Dr(fgm_,controlVariables_,variables_);

        for (int i = fgm_->Ncv; i < fgm_->Nvar; i++)
        {
            //Info << fgm_->varname[i] << endl;
            if (std::strcmp(fgm_->varname[i], "SOURCE_CV1") == 0){
                sourcePVCells[celli] = variables_[i];
                //Info << "sourcePVCells :" << sourcePVCells[celli] << endl;
            }
            if (std::strcmp(fgm_->varname[i], "PREFDIFF_CV1_H") == 0){
                prefDPV_HCells[celli] = variables_[i];
            }
            if (std::strcmp(fgm_->varname[i], "PREFDIFF_CV1_H2") == 0){
                prefDPV_H2Cells[celli] = variables_[i];
            }
            if (std::strcmp(fgm_->varname[i], "PREFDIFF_CV1_H2OPREFDIFF_CV1_T") == 0){ // BUG in variable read
                prefDPV_H2OCells[celli] = variables_[i];
                //Info << "prefDPV_H2OCells :" << prefDPV_H2OCells[celli] << endl;
            }
            if (std::strcmp(fgm_->varname[i], "PREFDIFF_CV1_T") == 0){
                prefDPV_TCells[celli] = variables_[i];
            }
            if (std::strcmp(fgm_->varname[i], "PREFDIFF_CV2_H") == 0){
                prefDh_HCells[celli] = variables_[i];
            }
            if (std::strcmp(fgm_->varname[i], "PREFDIFF_CV2_H2") == 0){
                prefDh_H2Cells[celli] = variables_[i];
            }
            if (std::strcmp(fgm_->varname[i], "PREFDIFF_CV2_H2OPREFDIFF_CV2_T") == 0){ // BUG in variable read
                prefDh_H2OCells[celli] = variables_[i];
                //Info << "prefDh_H2OCells :" << prefDh_H2OCells[celli] << endl;
            }
            if (std::strcmp(fgm_->varname[i], "PREFDIFF_CV2_T") == 0){
                prefDh_TCells[celli] = variables_[i];
            }
            if (std::strcmp(fgm_->varname[i], "PREFDIFF_CV3_H") == 0){
                prefDZ_HCells[celli] = variables_[i];
            }
            if (std::strcmp(fgm_->varname[i], "PREFDIFF_CV3_H2") == 0){
                prefDZ_H2Cells[celli] = variables_[i];
            }
            if (std::strcmp(fgm_->varname[i], "PREFDIFF_CV3_H2OPREFDIFF_CV3_T") == 0){ // BUG in variable read
                prefDZ_H2OCells[celli] = variables_[i];
                //Info << "prefDZ_H2OCells :" << prefDZ_H2OCells[celli] << endl;
            }
            if (std::strcmp(fgm_->varname[i], "PREFDIFF_CV3_T") == 0){
                prefDZ_TCells[celli] = variables_[i];
            }
            if (std::strcmp(fgm_->varname[i], "H") == 0){
                HCells[celli] = variables_[i];
            }
            if (std::strcmp(fgm_->varname[i], "H2") == 0){
                H2Cells[celli] = variables_[i];
            }
            if (std::strcmp(fgm_->varname[i], "H2O") == 0){
                H2OCells[celli] = variables_[i];
            }
            if (std::strcmp(fgm_->varname[i], "TEMPERATURE") == 0){
                TCells[celli] = variables_[i];
            }
            if (std::strcmp(fgm_->varname[i], "VISCOSITY") == 0){
                muCells[celli] = variables_[i];
            }
            if (std::strcmp(fgm_->varname[i], "CP") == 0){
                cp = variables_[i];
            }
            if (std::strcmp(fgm_->varname[i], "CONDUCTIVITY") == 0){
                lambda = variables_[i];
            }
            if (std::strcmp(fgm_->varname[i], "DENSITY") == 0){
                psiCells[celli] = variables_[i]/101325.0;
            }
        }

        // If cp and lambda are in the variable list, update alpha
        if (cp != 1.0e-10 && lambda != 1.0e-10){
            alphaCells[celli] = lambda/cp;
        }
        else {
            Info << "WARNING: Cp and/or lambda do not have a value, alpha does not have a value" << endl;
        }

        // Source term is strictly positive
        /*if (sourcePVCells[celli] < 0.0)
        {
            sourcePVCells[celli]   = 0.0;
        }
        if (PVsourcePVCells[celli] < 0.0)
        {
            PVsourcePVCells[celli]   = 0.0;
        }*/
    }

    forAll(T_.boundaryFieldRef(), patchi)
    {
        fvPatchScalarField& pZ           = this->Z_.boundaryFieldRef()[patchi];
        fvPatchScalarField& ph           = this->h_.boundaryFieldRef()[patchi];
        fvPatchScalarField& pPV          = this->PV_.boundaryFieldRef()[patchi];
        
        fvPatchScalarField& psourcePV    = this->sourcePV_.boundaryFieldRef()[patchi];
        fvPatchScalarField& pprefDPV_H   = this->prefDPV_H_.boundaryFieldRef()[patchi]; 
        fvPatchScalarField& pprefDPV_H2  = this->prefDPV_H2_.boundaryFieldRef()[patchi]; 
        fvPatchScalarField& pprefDPV_H2O = this->prefDPV_H2O_.boundaryFieldRef()[patchi]; 
        fvPatchScalarField& pprefDPV_T   = this->prefDPV_T_.boundaryFieldRef()[patchi]; 
        fvPatchScalarField& pprefDh_H    = this->prefDh_H_.boundaryFieldRef()[patchi]; 
        fvPatchScalarField& pprefDh_H2   = this->prefDh_H2_.boundaryFieldRef()[patchi]; 
        fvPatchScalarField& pprefDh_H2O  = this->prefDh_H2O_.boundaryFieldRef()[patchi]; 
        fvPatchScalarField& pprefDh_T    = this->prefDh_T_.boundaryFieldRef()[patchi]; 
        fvPatchScalarField& pprefDZ_H    = this->prefDZ_H_.boundaryFieldRef()[patchi]; 
        fvPatchScalarField& pprefDZ_H2   = this->prefDZ_H2_.boundaryFieldRef()[patchi]; 
        fvPatchScalarField& pprefDZ_H2O  = this->prefDZ_H2O_.boundaryFieldRef()[patchi]; 
        fvPatchScalarField& pprefDZ_T    = this->prefDZ_T_.boundaryFieldRef()[patchi]; 
        fvPatchScalarField& pH           = this->H_.boundaryFieldRef()[patchi];
        fvPatchScalarField& pH2          = this->H2_.boundaryFieldRef()[patchi];
        fvPatchScalarField& pH2O         = this->H2O_.boundaryFieldRef()[patchi];
        fvPatchScalarField& pT           = this->T_.boundaryFieldRef()[patchi];
        fvPatchScalarField& pmu          = this->mu_.boundaryFieldRef()[patchi];
        fvPatchScalarField& palpha       = this->alpha_.boundaryFieldRef()[patchi];
        fvPatchScalarField& ppsi         = this->psi_.boundaryFieldRef()[patchi];

        forAll(pT, facei)
        {
            //controlVariables_[0] = pZ[facei];
            //controlVariables_[1] = ph[facei];
            //controlVariables_[2] = pPV[facei];
            controlVariables_[0] = pPV[facei];
            controlVariables_[1] = ph[facei];
            controlVariables_[2] = pZ[facei];
 
            //lookupFGM_ND(fgm_,controlVariables_,variables_);
            lookupFGM_3Dr(fgm_,controlVariables_,variables_);
            
            for (int i = fgm_->Ncv; i < fgm_->Nvar; i++)
            {
                if (std::strcmp(fgm_->varname[i], "PREFDIFF_CV1_H") == 0){
                    pprefDPV_H[facei] = variables_[i];
                }
                if (std::strcmp(fgm_->varname[i], "PREFDIFF_CV1_H2") == 0){
                    pprefDPV_H2[facei] = variables_[i];
                }
                if (std::strcmp(fgm_->varname[i], "PREFDIFF_CV1_H2OPREFDIFF_CV1_T") == 0){ // BUG in variable read
                    pprefDPV_H2O[facei] = variables_[i];
                }
                if (std::strcmp(fgm_->varname[i], "PREFDIFF_CV1_T") == 0){
                    pprefDPV_T[facei] = variables_[i];
                }
                if (std::strcmp(fgm_->varname[i], "PREFDIFF_CV2_H") == 0){
                    pprefDh_H[facei] = variables_[i];
                }
                if (std::strcmp(fgm_->varname[i], "PREFDIFF_CV2_H2") == 0){
                    pprefDh_H2[facei] = variables_[i];
                }
                if (std::strcmp(fgm_->varname[i], "PREFDIFF_CV2_H2OPREFDIFF_CV2_T") == 0){ // BUG in variable read
                    pprefDh_H2O[facei] = variables_[i];
                }
                if (std::strcmp(fgm_->varname[i], "PREFDIFF_CV2_T") == 0){
                    pprefDh_T[facei] = variables_[i];
                }
                if (std::strcmp(fgm_->varname[i], "PREFDIFF_CV3_H") == 0){
                    pprefDZ_H[facei] = variables_[i];
                }
                if (std::strcmp(fgm_->varname[i], "PREFDIFF_CV3_H2") == 0){
                    pprefDZ_H2[facei] = variables_[i];
                }
                if (std::strcmp(fgm_->varname[i], "PREFDIFF_CV3_H2OPREFDIFF_CV3_T") == 0){ // BUG in variable read
                    pprefDZ_H2O[facei] = variables_[i];
                }
                if (std::strcmp(fgm_->varname[i], "PREFDIFF_CV3_T") == 0){
                    pprefDZ_T[facei] = variables_[i];
                }
                if (std::strcmp(fgm_->varname[i], "H") == 0){
                    pH[facei] = variables_[i];
                }
                if (std::strcmp(fgm_->varname[i], "H2") == 0){
                    pH2[facei] = variables_[i];
                }
                if (std::strcmp(fgm_->varname[i], "H2O") == 0){
                    pH2O[facei] = variables_[i];
                }
                if (std::strcmp(fgm_->varname[i], "TEMPERATURE") == 0){
                    if (!pT.fixesValue()){
                        pT[facei] = variables_[i];
                    }
                }
                if (std::strcmp(fgm_->varname[i], "VISCOSITY") == 0){
                    pmu[facei] = variables_[i];
                }
                if (std::strcmp(fgm_->varname[i], "CP") == 0){
                    cp = variables_[i];
                }
                if (std::strcmp(fgm_->varname[i], "CONDUCTIVITY") == 0){
                    lambda = variables_[i];
                }
                if (std::strcmp(fgm_->varname[i], "DENSITY") == 0){
                    ppsi[facei] = variables_[i]/101325.0;
                }
            }
            
            // No reaction at the boundary faces
            psourcePV[facei]   = 0.0;
            
            // If cp and lambda are in the variable list, update alpha
            if (cp != 0 && lambda != 0){
                palpha[facei] = lambda/cp;
            }
            else {
                Info << "WARNING: Cp and/or lambda do not have a value, alpha does not have a value" << endl;
            }
        }
    }
    
    Info << "Parameter update time: (" << clockTime_.timeIncrement() << " s)" << endl;
}


template<class ReactionThermo>
void Foam::combustionModels::baseFGM<ReactionThermo>::correct()
{
    update();
}


template<class ReactionThermo>
Foam::tmp<Foam::fvScalarMatrix>
Foam::combustionModels::baseFGM<ReactionThermo>::R
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
Foam::combustionModels::baseFGM<ReactionThermo>::Qdot() const
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
