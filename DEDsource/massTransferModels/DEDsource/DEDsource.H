/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2021 OpenCFD Ltd.
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
    Foam::meltingEvaporationModels::DEDsource

Description
    Mass transfer DEDsource model. Simple model driven by field value difference as:

    \f[
        \dot{m} = C \rho \alpha (T - T_{activate})/T_{activate}
    \f]

    where C is a model constant.

    if C > 0:
    \f[
        \dot{m} = C \rho \alpha (T - T_{activate})/T_{activate}
    \f]
        for \f[ T > T_{activate} \f]

    and

        \f[ mDot = 0.0 \f]  for \f[  T < T_{activate} \f]


    if C < 0:
    \f[
        \dot{m} = -C \rho \alpha (T_{activate} - T)/T_{activate}
    \f]
        for \f[ T < T_{activate} \f]

    and
        \f[ \dot{m} = 0.0 \f] for \f[ T > T_{activate} \f]

    Based on the reference:
    -# W. H. DEDsource. "A Pressure Iteration Scheme for Two-Phase Modeling".
    Technical Report LA-UR 79-975. Los Alamos Scientific Laboratory,
    Los Alamos, New Mexico. 1979.

Usage
    Example usage:
    \verbatim
        massTransferModel
        (
            (solid to liquid)
            {
                type            DEDsource;
                C               40;
                Tactivate       302.78;
            }
        );
    \endverbatim

    Where:

    \table
        Property    | Description            | Required    | Default value
        Tactivate   | Activation temperature | yes
        C           | Model constant         | yes
        includeVolChange    | Volumen change  | no          | yes
        species     | Specie name on the other phase | no   | none
    \endtable

SourceFiles
    DEDsource.C

\*---------------------------------------------------------------------------*/

#ifndef meltingEvaporationModels_DEDsource_H
#define meltingEvaporationModels_DEDsource_H

#include "InterfaceCompositionModel.H"

// * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * * * * *//

namespace Foam
{
namespace meltingEvaporationModels
{

/*---------------------------------------------------------------------------*\
                             Class DEDsource Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo, class OtherThermo>
class DEDsource
:
    public InterfaceCompositionModel<Thermo, OtherThermo>
{
    // Private Data

		//-Radius of powder nozzle[m]
		dimensionedScalar Rnozzle_;
		
		//-Working height [m]
		dimensionedScalar Hworking_;
		
		//-Injection angle []
		dimensionedScalar angleInj_;
		
		//-Divergence angle []
		dimensionedScalar angleDiv_;
		
		//-Start of injection [s]
		dimensionedScalar SOI_;
		
		//-Duration [s]
		dimensionedScalar duration_;
		
		//-Starting position [m]
		dimensionedVector position0_;
		
		//-Processing velocity [m/s]
		dimensionedVector  vel_;
		
		//-Powder flowrate [kg/s]
		dimensionedScalar mdot_;
		
		//-Powder distribution field [kg/m3/s]
		volScalarField powderDistribution_;
		
        //- Condensation coefficient [1/s]
        dimensionedScalar C_;

        //- Phase transition temperature
        const dimensionedScalar Tactivate_;

        //- Phase minimum value for activation
        scalar alphaMin_;


public:

    //- Runtime type information
    TypeName("DEDsource");


    // Constructors

        //- Construct from components
        DEDsource
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~DEDsource() = default;


    // Member Functions

        //- Explicit total mass transfer coefficient
        virtual tmp<volScalarField> Kexp
        (
            const volScalarField& field
        );

        //- Implicit mass transfer coefficient
        virtual tmp<volScalarField> KSp
        (
            label modelVariable,
            const volScalarField& field
        );

        //- Explicit mass transfer coefficient
        virtual tmp<volScalarField> KSu
        (
            label modelVariable,
            const volScalarField& field
        );

        //- Return T transition between phases
        virtual const dimensionedScalar& Tactivate() const noexcept
        {
            return Tactivate_;
        }

        //- Add/subtract alpha*div(U) as a source term
        //- for alpha, substituting div(U) = mDot(1/rho1 - 1/rho2)
        virtual bool includeDivU() const noexcept
        {
            return true;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace meltingEvaporationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "DEDsource.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
