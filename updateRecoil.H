psat = pstd*Foam::exp(1e-3*mw*L_vapor/R*(1/Tv - 1/T));
// forAll(C, i)
// {
	// if(T[i] <= T_V.value())
	// {
		// psat[i] = 0.0;
	// }
// }
psat.correctBoundaryConditions();

recoil = psat*fvc::grad(alpha1)*factor;
recoil.correctBoundaryConditions();
