/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2022 OpenCFD Ltd.
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

#include "DEDsource.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo, class OtherThermo>
Foam::meltingEvaporationModels::DEDsource<Thermo, OtherThermo>::DEDsource
(
    const dictionary& dict,
    const phasePair& pair
)
:
    InterfaceCompositionModel<Thermo, OtherThermo>(dict, pair),
	Rnozzle_("Rnozzle", dimLength, dict),
	Hworking_("Hworking", dimLength, dict),
	angleInj_("angleInj", dimless, dict),
	angleDiv_("angleDiv", dimless, dict),
	SOI_("SOI", dimTime, dict),
	duration_("duration", dimTime, dict),
	position0_("StartingPosition", dimLength, dict),
	vel_("ProcessVelocity", dimVelocity, dict),
	mdot_("MassFlowrate", dimensionSet(1, 0, -1, 0, 0, 0, 0), dict),
	powderDistribution_
	(
		IOobject
        (
            "powderDistribution",
            this->mesh_.time().timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar(dimDensity/dimTime, Zero)
	),
    C_("C", inv(dimTime), dict),
    Tactivate_("Tactivate", dimTemperature, dict),
    alphaMin_(dict.getOrDefault<scalar>("alphaMin", 0))
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Thermo, class OtherThermo>
Foam::tmp<Foam::volScalarField>
Foam::meltingEvaporationModels::DEDsource<Thermo, OtherThermo>::Kexp
(
    const volScalarField& refValue
)
{
	const fvMesh& mesh = this->mesh_;
	const volVectorField& C = mesh.C();
	const volScalarField& rho1 = this->pair().to().rho();
	const volScalarField& rho2 = this->pair().from().rho();
	const volScalarField& alpha1 = this->pair().to();
	const volScalarField& alpha2 = this->pair().from();
	const volScalarField rho = alpha1*rho1 + alpha2*rho2;

	auto tdelta = tmp<volScalarField>::New
	(
		IOobject
		(
			"tdelta  ",
			mesh.time().timeName(),
			mesh
		),
		mesh,
		dimensionedScalar(inv(dimLength), 0.0)
	);
	volScalarField& delta  = tdelta.ref();
	
	delta = 2*rho/(rho1 + rho2)*mag(Foam::fvc::grad(alpha1));
	
	if(mesh.time().value() >=  SOI_.value() && mesh.time().value() <=  (SOI_ + duration_).value())
	{
		auto tp0 = tmp<volVectorField>::New
		(
			IOobject
			(
				"tp0 ",
				mesh.time().timeName(),
				mesh
			),
			mesh,
			dimensionedVector(dimLength, vector(0.0, 0.0, 0.0))
		);
		volVectorField& p0  = tp0.ref();
		
		p0 = C - position0_;
		
		dimensionedScalar injAngle = angleInj_*Foam::constant::mathematical::pi / 180.0;
		dimensionedScalar divAngle = angleDiv_*Foam::constant::mathematical::pi / 180.0;
		
		{
			auto tnozzlePos1 = tmp<volVectorField>::New
			(
				IOobject
				(
					"tnozzlePos1",
					mesh.time().timeName(),
					mesh
				),
				mesh,
				dimensionedVector(dimLength, vector(0.0, 0.0, 0.0))
			);
			volVectorField& nozzlePos1 = tnozzlePos1.ref();
		
			nozzlePos1.component(0) = (p0.component(0) - Hworking_/Foam::tan(injAngle))*Foam::sin(injAngle) - (p0.component(1) - Hworking_)*Foam::cos(injAngle);
			nozzlePos1.component(1) = (p0.component(0) - Hworking_/Foam::tan(injAngle))*Foam::cos(injAngle) + (p0.component(1) - Hworking_)*Foam::sin(injAngle);
			nozzlePos1.component(2) = p0.component(2);
			
			auto tnozzlePos2 = tmp<volVectorField>::New
			(
				IOobject
				(
					"tnozzlePos2",
					mesh.time().timeName(),
					mesh
				),
				mesh,
				dimensionedVector(dimLength, vector(0.0, 0.0, 0.0))
			);
			volVectorField& nozzlePos2 = tnozzlePos2.ref();
		
			nozzlePos2.component(0) = p0.component(0) ;
			nozzlePos2.component(1) = (p0.component(2) - Hworking_/Foam::tan(injAngle))*Foam::cos(injAngle) + (p0.component(1) - Hworking_)*Foam::sin(injAngle);
			nozzlePos2.component(2) = (p0.component(2) - Hworking_/Foam::tan(injAngle))*Foam::sin(injAngle) - (p0.component(1) - Hworking_)*Foam::cos(injAngle);
			
		
			auto tnozzlePos3 = tmp<volVectorField>::New
			(
				IOobject
				(
					"tnozzlePos3",
					mesh.time().timeName(),
					mesh
				),
				mesh,
				dimensionedVector(dimLength, vector(0.0, 0.0, 0.0))
			);
			volVectorField& nozzlePos3 = tnozzlePos3.ref();
		
			nozzlePos3.component(0) = (p0.component(0) + Hworking_/Foam::tan(injAngle))*Foam::sin(injAngle) + (p0.component(1) - Hworking_)*Foam::cos(injAngle);
			nozzlePos3.component(1) = - ((p0.component(0) + Hworking_/Foam::tan(injAngle))*Foam::cos(injAngle) - (p0.component(1) - Hworking_)*Foam::sin(injAngle));
			nozzlePos3.component(2) = p0.component(2);
			
			auto tnozzlePos4 = tmp<volVectorField>::New
			(
				IOobject
				(
					"tnozzlePos4",
					mesh.time().timeName(),
					mesh
				),
				mesh,
				dimensionedVector(dimLength, vector(0.0, 0.0, 0.0))
			);
			volVectorField& nozzlePos4 = tnozzlePos4.ref();
		
			nozzlePos4.component(0) = p0.component(0) ;
			nozzlePos4.component(1) = -((p0.component(2)+ Hworking_/Foam::tan(injAngle))*Foam::cos(injAngle) - (p0.component(1) - Hworking_)*Foam::sin(injAngle));
			nozzlePos4.component(2) = (p0.component(2) + Hworking_/Foam::tan(injAngle))*Foam::sin(injAngle) + (p0.component(1) - Hworking_)*Foam::cos(injAngle);

		
			volScalarField  r1 = Rnozzle_ - nozzlePos1.component(1)*Foam::tan(divAngle);
			volScalarField  r2 = Rnozzle_ - nozzlePos2.component(1)*Foam::tan(divAngle);
			volScalarField  r3 = Rnozzle_ - nozzlePos3.component(1)*Foam::tan(divAngle);
			volScalarField  r4 = Rnozzle_ - nozzlePos4.component(1)*Foam::tan(divAngle);
		
			volScalarField  A1 = Foam::constant::mathematical::pi*Foam::sqr(r1);
			volScalarField  A2 = Foam::constant::mathematical::pi*Foam::sqr(r2);
			volScalarField  A3 = Foam::constant::mathematical::pi*Foam::sqr(r3);
			volScalarField  A4 = Foam::constant::mathematical::pi*Foam::sqr(r4);
		

			auto tS1 = tmp<volScalarField>::New
			(
				IOobject
				(
					"tS1 ",
					mesh.time().timeName(),
					mesh
				),
				mesh,
				dimensionedScalar(dimDensity/dimTime*dimLength, Zero)
			);
			volScalarField& S1  = tS1.ref();
			
			S1 = 0.5*mdot_/A1*Foam::exp(-2*(Foam::sqr(nozzlePos1.component(0)) + Foam::sqr(nozzlePos1.component(2))) / Foam::sqr(r1));
			
			auto tS2 = tmp<volScalarField>::New
			(
				IOobject
				(
					"tS2 ",
					mesh.time().timeName(),
					mesh
				),
				mesh,
				dimensionedScalar(dimDensity/dimTime*dimLength, Zero)
			);
			volScalarField& S2  = tS2.ref();
			
			S2 = 0.5*mdot_/A2*Foam::exp(-2*(Foam::sqr(nozzlePos2.component(0)) + Foam::sqr(nozzlePos2.component(2))) / Foam::sqr(r2));
			
			auto tS3 = tmp<volScalarField>::New
			(
				IOobject
				(
					"tS3 ",
					mesh.time().timeName(),
					mesh
				),
				mesh,
				dimensionedScalar(dimDensity/dimTime*dimLength, Zero)
			);
			volScalarField& S3  = tS3.ref();
			
			S3 = 0.5*mdot_/A3*Foam::exp(-2*(Foam::sqr(nozzlePos3.component(0)) + Foam::sqr(nozzlePos3.component(2))) / Foam::sqr(r3));
			
			auto tS4 = tmp<volScalarField>::New
			(
				IOobject
				(
					"tS4 ",
					mesh.time().timeName(),
					mesh
				),
				mesh,
				dimensionedScalar(dimDensity/dimTime*dimLength, Zero)
			);
			volScalarField& S4  = tS4.ref();
			
			S4 = 0.5*mdot_/A4*Foam::exp(-2*(Foam::sqr(nozzlePos4.component(0)) + Foam::sqr(nozzlePos4.component(2))) / Foam::sqr(r4));
			
			powderDistribution_ = pos0(refValue - Tactivate_)*(S1 + S2 + S3 + S4)*delta;
		}
	}
	else
	{
		powderDistribution_  = dimensionedScalar(dimDensity/dimTime,  Zero);
	}
	
    return tmp<volScalarField>::New(powderDistribution_);
}


template<class Thermo, class OtherThermo>
Foam::tmp<Foam::volScalarField>
Foam::meltingEvaporationModels::DEDsource<Thermo, OtherThermo>::KSp
(
    label variable,
    const volScalarField& refValue
)
{
	return tmp<volScalarField> ();
}


template<class Thermo, class OtherThermo>
Foam::tmp<Foam::volScalarField>
Foam::meltingEvaporationModels::DEDsource<Thermo, OtherThermo>::KSu
(
    label variable,
    const volScalarField& refValue
)
{
	const volScalarField& rhoTo = this->pair().to().rho();
	const volScalarField& rhoFrom = this->pair().from().rho();

	const volScalarField& alphaTo = this->pair().to();
	const volScalarField& alphaFrom = this->pair().from();
	
	const volScalarField rho = alphaTo*rhoTo + alphaFrom*rhoFrom;
	
	if(interfaceCompositionModel::P == variable)
	{
		volScalarField coeff(powderDistribution_*pos(alphaTo - 0.9)*(rhoFrom/(rhoFrom - rhoTo)));
		
		return tmp<volScalarField>::New(coeff);
	}
	else	if(interfaceCompositionModel::alpha == variable )
	{
		volScalarField coeff(powderDistribution_*rhoFrom/rho*pos(refValue - 0.9));		
		
		return tmp<volScalarField>::New(coeff);
	}
	else
		return tmp<volScalarField> ();
}


// ************************************************************************* //
