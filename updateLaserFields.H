//catch the center point of each mesh
Info << "++Updating laser motion." << endl;


if(runTime.value() >= startTime.value() && runTime.value() <= (startTime + duration).value()) //start the laser until the flowtime gets the laser begin time
{
	volScalarField r_local = Foam::sqrt(
									 Foam::sqr(C.component(0) - X0.component(0) - Ulaser.component(0)*(runTime.time() - startTime))
									+Foam::sqr(C.component(2) - X0.component(2) - Ulaser.component(2)*(runTime.time() - startTime))
									); //calculate the laser
	Q_laser = eta*scalar(2.0)*Plaser/(Foam::constant::mathematical::pi * Foam::sqr(Rlaser)) * Foam::exp(-2 * Foam::sqr(r_local/Rlaser)) * mag(fvc::grad(alpha1))*factor;
}
else
	Q_laser = dimensionedScalar(dimensionSet(1, -1, -3, 0, 0, 0, 0), Zero);

Q_laser.correctBoundaryConditions();
