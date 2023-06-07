Info << "Reading fusion heat properties from file.\n" << endl;
IOdictionary fusionProperties//OpenFOAM dinctionary
(
	IOobject
	(
		"fusionProperties",
		runTime.constant(),
		mesh,
		IOobject::MUST_READ,
		IOobject::NO_WRITE
	)
);

//Solidus Temp. for alpha1
dimensionedScalar Ts
(
	"Ts",
	dimensionSet(0, 0, 0, 1, 0, 0, 0),
	fusionProperties,
	"Ts"
);


//Liquidus Temp. for alpha1
dimensionedScalar Tl
(
	"Tl",
	dimensionSet(0, 0, 0, 1, 0, 0, 0),
	fusionProperties,
	"Tl"
);

//Vaporization Temp.
dimensionedScalar Tv
(
	"Tv",
	dimensionSet(0, 0, 0, 1, 0, 0, 0),
	fusionProperties,
	"Tv"
);

//Fusion heat for alpha1
dimensionedScalar L_fusion
(
	"L_fusion",
	dimensionSet(0, 2, -2, 0, 0, 0, 0),
	fusionProperties,
	"L_fusion"
);

dimensionedScalar L_vapor
(
	"L_vapor",
	dimensionSet(0, 2, -2, 0, 0, 0, 0),
	fusionProperties,
	"L_vapor"
);

dimensionedScalar mw
(
	"mw",
	dimensionSet(1, 0, 0, 0, -1, 0, 0),
	fusionProperties,
	"mw"
);

Info << "++Solidus Temperature1:" << Ts << endl;
Info << "++Liquidus Temperature1:" << Tl << endl;
Info << "++L_fusion1:" << L_fusion << endl;
Info << "++L_vapor:" << L_vapor << endl;
