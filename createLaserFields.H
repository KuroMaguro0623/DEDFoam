Info << "Reading laser properties from file.\n" << endl;
IOdictionary beamProperties//OpenFOAM dinctionary
(
	IOobject
	(
		"beamProperties",
		runTime.constant(),
		mesh,
		IOobject::MUST_READ,
		IOobject::NO_WRITE
	)
);

//laser speed
dimensionedVector Ulaser
(
	"Ulaser",
	dimensionSet(0, 1, -1, 0, 0, 0, 0),
	beamProperties,
	"Ulaser"
);


//laser initial position
dimensionedVector X0
(
	"X0",
	dimensionSet(0, 1, 0, 0, 0, 0, 0),
	beamProperties,
	"X0"
);


//laser radius
dimensionedScalar Rlaser
(
	"Rlaser",
	dimensionSet(0, 1, 0, 0, 0, 0, 0),
	beamProperties,
	"Rlaser"
);

//laser power
dimensionedScalar Plaser
(
	"Plaser",
	dimensionSet(1, 2, -3, 0, 0, 0, 0),
	beamProperties,
	"Plaser"
);

//absorbtivity
dimensionedScalar eta
(
	"eta",
	dimensionSet(0, 0, 0, 0, 0, 0, 0),
	beamProperties,
	"eta"
);

//startTime
dimensionedScalar startTime
(
	"startTime",
	dimensionSet(0, 0, 1, 0, 0, 0, 0),
	beamProperties,
	"startTime"
);

//duration
dimensionedScalar duration
(
	"duration",
	dimensionSet(0, 0, 1, 0, 0, 0, 0),
	beamProperties,
	"duration"
);

Info << "++Ulaser:" << Ulaser << endl;
Info << "++Initial position:" << X0 << endl;
Info << "++Laser radius:" << Rlaser << endl;
Info << "++Laser power:" << Plaser << endl;
Info << "++Laser absorptivity:" << eta << endl;
Info << "++Laser time to begin:" << startTime << endl;
Info << "++Laser duration:" << duration << endl;

volScalarField Q_laser//laser source term
(
	IOobject
	(
		"Q_laser",
		runTime.timeName(),
		mesh,
		IOobject::READ_IF_PRESENT,
		IOobject::AUTO_WRITE
	),
	T.mesh(),
	dimensionedScalar
	(
		"Q_laser",
		dimensionSet(1, -1, -3, 0, 0, 0, 0),
		scalar(0.0)
	
	),
	zeroGradientFvPatchScalarField::typeName
);
