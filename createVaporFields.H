Info << "++Vapor parameters are created." << endl;

const dimensionedScalar pstd("pstd", dimPressure, Foam::constant::thermodynamic::Pstd);
const dimensionedScalar R = Foam::constant::physicoChemical::R;


volScalarField psat
(
    IOobject
    (
        "psat",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
	p.mesh(),
	dimensionedScalar
	(
		"psat",
		dimPressure,
		scalar(0.0)
	),
	zeroGradientFvPatchScalarField::typeName
);

volVectorField recoil
(
    IOobject
    (
        "recoil",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
	psat*fvc::grad(alpha1)
);

volScalarField Q_vapor
(
	IOobject
	(
		"Q_vapor",
		runTime.timeName(),
		mesh,
		IOobject::READ_IF_PRESENT,
		IOobject::AUTO_WRITE
	),
	T.mesh(),
	dimensionedScalar
	(
		"Q_vapor",
		dimensionSet(1, -1, -3, 0, 0, 0, 0),
		scalar(0.0)
	),
	zeroGradientFvPatchScalarField::typeName
);

volScalarField m_vap
(
	IOobject
	(
		"m_vap",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	),
	T.mesh(),
	dimensionedScalar
	(
		"m_vap",
		dimensionSet(1, -2, -1, 0, 0, 0, 0),
		scalar(0.0)
	),
	zeroGradientFvPatchScalarField::typeName
);
