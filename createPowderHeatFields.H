Info << "Reading powder heat properties from file.\n" << endl;
IOdictionary powderHeatProperties//OpenFOAM dinctionary
(
	IOobject
	(
		"powderHeatProperties",
		runTime.constant(),
		mesh,
		IOobject::MUST_READ,
		IOobject::NO_WRITE
	)
);

//powderHeat
dimensionedScalar h_powder
(
	"h_powder",
	dimensionSet(0, 2,-2, 0, 0, 0, 0),
	powderHeatProperties,
	"h_powder"
);

volScalarField Q_powder//Powder source term
(
	IOobject
	(
		"Q_powder",
		runTime.timeName(),
		mesh,
		IOobject::READ_IF_PRESENT,
		IOobject::AUTO_WRITE
	),
	mesh,
	dimensionedScalar
	(
		"Q_powder",
		dimensionSet(1, -1, -3, 0, 0, 0, 0),
		scalar(0.0)
	)
);
