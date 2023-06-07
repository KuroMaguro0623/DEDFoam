volScalarField beta
(
    IOobject
    (
        "beta",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    T.mesh(),
	dimensionedScalar
	(
		"beta",
		dimless,
		scalar(0.0)
	),
	zeroGradientFvPatchScalarField::typeName
);

volScalarField A_c
(
    IOobject
    (
        "A_c",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    T.mesh(),
	dimensionedScalar
	(
		"A_c",
		dimensionSet(1, -3, -1, 0, 0, 0, 0),
		scalar(0.0)
	),
	zeroGradientFvPatchScalarField::typeName
);

dimensionedScalar A_mush
(
	"A_mush",
	dimensionSet(1, -3, -1, 0, 0, 0, 0),
	scalar(1e12)
);

dimensionedScalar bb
(
	"bb",
	dimless,
	scalar(1e-3)
);

