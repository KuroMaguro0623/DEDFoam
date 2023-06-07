Info << "Reading properties for Marangoni effects.\n" << endl;

volVectorField Marangoni
(
	IOobject
    (
		"Marangoni",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
    ),
	mesh,
	dimensionedVector
	(
		"Marangoni",
		dimensionSet(1, -2, -2, 0, 0, 0, 0),
		vector(0.0, 0.0, 0.0)
	)
);

IOdictionary MarangoniProperties
(
	IOobject
	(
		"MarangoniProperties",
		runTime.constant(),
		mesh,
		IOobject::MUST_READ,
		IOobject::NO_WRITE
	)
);

dimensionedScalar dSigma_dT
(
	"dSigma_dT",
	dimForce/dimLength/dimTemperature,
	MarangoniProperties,
	"dSigma_dT"
);

