/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

Application
    postLES

Description
    Calculates and writes some mean fields from an LES calculation. The
    following fields are created:
        - resLES: resolutness of an LES
        - TKEMean: mean turbulent kinectic energy
        - TKEMeanProd: production term of the mean turbulent kinectic energy
		- turbDiffusionMean: turbulent diffusion term of the mean turbulent kinectic energy
		- SGSDiffusionMean: subgrid scale diffusion term of the mean turbulent kinectic energy
		- viscDiffusionMean: viscous diffusion term of the mean turbulent kinectic energy
    Fields UMean, UPrime2Mean, turbDiffMean, SGSDiffMean, and kMean must exist.

\*---------------------------------------------------------------------------*/

// #include "calc.H"
#include "fvc.H"
#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();
    timeSelector::addOptions();

    #include "setRootCase.H"
    #include "createTime.H"

    // Get times list
    instantList timeDirs = timeSelector::select0(runTime, args);

    #include "createMesh.H"

    // new
    #include "createControl.H"
    #include "createTimeControls.H"
	
    forAll(timeDirs, timeI)
    {   
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;	
		
		IOobject UPrime2MeanHeader
		(
			"UPrime2Mean",
			runTime.timeName(),
			mesh,
			IOobject::MUST_READ
		);

		IOobject UMeanHeader
		(
			"UMean",
			runTime.timeName(),
			mesh,
			IOobject::MUST_READ
		);
			  
		IOobject convectionMeanHeader
		(
			"convectionMean",
			runTime.timeName(),
			mesh,
			IOobject::MUST_READ
		);         
		
		IOobject transportTurbulenceFieldMeanHeader
		(
			"transportTurbulenceFieldMean",
			runTime.timeName(),
			mesh,
			IOobject::MUST_READ
		);     
		
		IOobject transportVelocityPressureFieldMeanHeader
		(
			"transportVelocityPressureFieldMean",
			runTime.timeName(),
			mesh,
			IOobject::MUST_READ
		); 
		
		IOobject viscousDissipationMeanHeader
		(
			"viscousDissipationMean",
			runTime.timeName(),
			mesh,
			IOobject::MUST_READ
		); 	
	
		if
		(
			(UPrime2MeanHeader.typeHeaderOk<volSymmTensorField>(true))
		 && (UMeanHeader.typeHeaderOk<volVectorField>(true))
		 && (convectionMeanHeader.typeHeaderOk<volScalarField>(true))
		 && (transportTurbulenceFieldMeanHeader.typeHeaderOk<volVectorField>(true))     
		 && (transportVelocityPressureFieldMeanHeader.typeHeaderOk<volVectorField>(true))	 
		 && (viscousDissipationMeanHeader.typeHeaderOk<volScalarField>(true))	 
		)		
		{
			
			Info<< "    Reading average field UPrime2Mean" << endl;
			const volSymmTensorField UPrime2Mean(UPrime2MeanHeader, mesh);

			Info<< "    Reading average field UMean" << endl;
			const volVectorField UMean(UMeanHeader, mesh);                

			Info<< "    Reading average field convectionMean" << endl;
			const volScalarField convectionMean(convectionMeanHeader, mesh);		
			
			Info<< "    Reading average field transportTurbulenceFieldMean" << endl;
			const volVectorField transportTurbulenceFieldMean(transportTurbulenceFieldMeanHeader, mesh);                                  
			  
			Info<< "    Reading average field transportVelocityPressureFieldMean" << endl;
			const volVectorField transportVelocityPressureFieldMean(transportVelocityPressureFieldMeanHeader, mesh);

			Info<< "    Reading average field transportVelocityPressureMean" << endl;
			const volScalarField viscousDissipationMean(viscousDissipationMeanHeader, mesh);		
			
			volVectorField U
			(
				IOobject
				(
					"U",
					runTime.timeName(),
					mesh,
					IOobject::MUST_READ,
					IOobject::NO_WRITE
				),
				mesh
			);                 


	///////////////////////////////////////////////////////////////////////        

			#include "createPhi.H"

			singlePhaseTransportModel laminarTransport(U, phi);

			autoPtr<incompressible::turbulenceModel> turbulence
			(
				incompressible::turbulenceModel::New(U, phi, laminarTransport)
			);
		
			volScalarField kMean
			(
				IOobject
				(
				"kMean",
				runTime.timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::AUTO_WRITE
				),
				0.5*tr(UPrime2Mean)
			);	
		
			volScalarField productionMean
			(
				IOobject
				(
				"productionMean",
				runTime.timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::AUTO_WRITE
				),
				-UPrime2Mean && fvc::grad(UMean)
			);
			
			volScalarField transportViscousDiffusionMean
			(
				IOobject
				(
				"transportViscousDiffusionMean",
				runTime.timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::AUTO_WRITE
				),
				turbulence->nuEff()*fvc::laplacian(kMean)
			);		
	 
			volScalarField transportVelocityPressureMean
			(
				IOobject
				(
				"transportVelocityPressureMean",
				runTime.timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::AUTO_WRITE
				),
				-fvc::div(transportVelocityPressureFieldMean)
			); 
			
			volScalarField transportTurbulenceMean
			(
				IOobject
				(
				"transportTurbulenceMean",
				runTime.timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::AUTO_WRITE
				),
				-fvc::div(transportTurbulenceFieldMean)
			);	
			
			productionMean.write();
			transportViscousDiffusionMean.write();
			transportVelocityPressureMean.write(); 
			transportTurbulenceMean.write();            
				   
		}
		else
		{
			Info << "    No UPrime2Mean and/or kMean and/or UMean." << endl;
		}
    }
    Info<< "End" << endl;
}


// ************************************************************************* //
