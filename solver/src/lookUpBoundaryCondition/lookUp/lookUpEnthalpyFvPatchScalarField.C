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
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "lookUpEnthalpyFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H" 
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::lookUpEnthalpyFvPatchScalarField::lookUpEnthalpyFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    fgm_(nullptr),
    iT(-1),
    icp(-1),
    h0_(p.size(), Zero)
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

    fgm_ = FGMTable::getFGMInstance("constant/lookUp/database.fgm");

    // Initialize arrays with sizes obtained from the FGM table
    controlVariables_ = new float[fgm_->Ncv];
    variables_        = new float[fgm_->Nvar];

    for (int i = fgm_->Ncv; i < fgm_->Nvar; i++)
    {
        if (std::strcmp(fgm_->varname[i], "TEMPERATURE") == 0){
            iT = i;
        }
        if (std::strcmp(fgm_->varname[i], "CP") == 0){
            icp = i;
        }
    }
}


Foam::lookUpEnthalpyFvPatchScalarField::lookUpEnthalpyFvPatchScalarField
(
    const lookUpEnthalpyFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    fgm_(nullptr),
    iT(-1),
    icp(-1),
    h0_(ptf.h0_, mapper)
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

    fgm_ = FGMTable::getFGMInstance("constant/lookUp/database.fgm");

    // Initialize arrays with sizes obtained from the FGM table
    controlVariables_ = new float[fgm_->Ncv];
    variables_        = new float[fgm_->Nvar];

    for (int i = fgm_->Ncv; i < fgm_->Nvar; i++)
    {
        if (std::strcmp(fgm_->varname[i], "TEMPERATURE") == 0){
            iT = i;
        }
        if (std::strcmp(fgm_->varname[i], "CP") == 0){
            icp = i;
        }
    }
}

Foam::lookUpEnthalpyFvPatchScalarField::lookUpEnthalpyFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict, false),
    fgm_(nullptr),
    iT(-1),
    icp(-1),
    h0_("h0", dict, p.size())
{
    if (dict.found("value"))
    {
        fvPatchField<scalar>::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        fvPatchField<scalar>::operator=(h0_);
    }
    
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

    fgm_ = FGMTable::getFGMInstance("constant/lookUp/database.fgm");
    		
    // Initialize arrays with sizes obtained from the FGM table
    controlVariables_ = new float[fgm_->Ncv];
    variables_        = new float[fgm_->Nvar];

    for (int i = fgm_->Ncv; i < fgm_->Nvar; i++)
    {
        if (std::strcmp(fgm_->varname[i], "TEMPERATURE") == 0){
            iT = i;
        }
        if (std::strcmp(fgm_->varname[i], "CP") == 0){
            icp = i;
        }
    }
}


Foam::lookUpEnthalpyFvPatchScalarField::lookUpEnthalpyFvPatchScalarField
(
    const lookUpEnthalpyFvPatchScalarField& tppsf
)
:
    fixedValueFvPatchScalarField(tppsf),
    fgm_(nullptr),
    iT(-1),
    icp(-1),
    h0_(tppsf.h0_)
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

    fgm_ = FGMTable::getFGMInstance("constant/lookUp/database.fgm");
    	
    // Initialize arrays with sizes obtained from the FGM table
    controlVariables_ = new float[fgm_->Ncv];
    variables_        = new float[fgm_->Nvar];

    for (int i = fgm_->Ncv; i < fgm_->Nvar; i++)
    {
        if (std::strcmp(fgm_->varname[i], "TEMPERATURE") == 0){
            iT = i;
        }
        if (std::strcmp(fgm_->varname[i], "CP") == 0){
            icp = i;
        }
    }
}

Foam::lookUpEnthalpyFvPatchScalarField::lookUpEnthalpyFvPatchScalarField
(
    const lookUpEnthalpyFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF),
    fgm_(nullptr),
    iT(-1),
    icp(-1),
    h0_(tppsf.h0_)
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

    fgm_ = FGMTable::getFGMInstance("constant/lookUp/database.fgm");
    		
    // Initialize arrays with sizes obtained from the FGM table
    controlVariables_ = new float[fgm_->Ncv];
    variables_        = new float[fgm_->Nvar];

    for (int i = fgm_->Ncv; i < fgm_->Nvar; i++)
    {
        if (std::strcmp(fgm_->varname[i], "TEMPERATURE") == 0){
            iT = i;
        }
        if (std::strcmp(fgm_->varname[i], "CP") == 0){
            icp = i;
        }
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::lookUpEnthalpyFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchScalarField::autoMap(m);
    h0_.autoMap(m);
}


void Foam::lookUpEnthalpyFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchScalarField::rmap(ptf, addr);

    const lookUpEnthalpyFvPatchScalarField& tiptf =
        refCast<const lookUpEnthalpyFvPatchScalarField>(ptf);

    h0_.rmap(tiptf.h0_, addr);
}


void Foam::lookUpEnthalpyFvPatchScalarField::updateCoeffs()
{

    if (updated())
    {
        return;
    }
    
    scalar T0, dha, cp, T1, Tresidual = 1.0;
    
    int iter;

    const fvPatchField<scalar>& pT =
        patch().lookupPatchField<volScalarField, scalar>("T");
			
    forAll(h0_, facei)
    {
        
        // Set control variables
	    controlVariables_[0] = *const_cast<scalar*>(&(patch().lookupPatchField<volScalarField, scalar>("Yc")[facei]));
	    controlVariables_[1] = *const_cast<scalar*>(&(patch().lookupPatchField<volScalarField, scalar>("h")[facei]));
		
        lookupFGM_2D(fgm_,controlVariables_,variables_);
		
        T0 = pT[facei];
        T1 = variables_[iT];
        cp = variables_[icp];

        int lookupStatus = EXIT_SUCCESS;
        iter = 0;

        while ((iter < 1000) &&
            (lookupStatus == EXIT_SUCCESS) &&
            (Foam::mag(T0 - T1) > Tresidual))
        {
            dha = 0.1*(T0 - T1) * cp;

            controlVariables_[1] += dha;

            if (lookupFGM_2D(fgm_,controlVariables_,variables_) == EXIT_SUCCESS)
            {
                T1 = variables_[iT];
                cp = variables_[icp];
            }
            iter++;
        }

        if (iter == 1000)
        {
            WarningInFunction
                << "WARNING: Iteration limit reached in lookupEnthalpyFvPatchScalarField::updateCoeffs()"
                << " for face " << facei
                << " with control variables: " << controlVariables_[0] << " "
                << controlVariables_[1] << endl;
        }
        else if (lookupStatus != EXIT_SUCCESS)
        {
            WarningInFunction
                << "WARNING: Lookup failed in lookupEnthalpyFvPatchScalarField::updateCoeffs()"
                << " for face " << facei
                << " with control variables: " << controlVariables_[0] << " " << controlVariables_[1]
                << " Temperature residual: " << Foam::mag(T0 - T1)
                << endl;
        }
	
        h0_[facei] = controlVariables_[1];
    }

    operator==
    (
        h0_
    );

    fixedValueFvPatchScalarField::updateCoeffs();	
}


void Foam::lookUpEnthalpyFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    h0_.writeEntry("h0", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        lookUpEnthalpyFvPatchScalarField
    );
}

// ************************************************************************* //
