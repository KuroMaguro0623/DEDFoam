Info << "++Updating fusion heat." << endl;

volScalarField cp_fusion = alpha1*L_fusion/(Tl - Ts);
forAll(C, i)
{
	if(T[i] >= Tl.value() || T[i] <= Ts.value())
	{
		cp_fusion[i] = 0.0;
	}
}
cp_fusion.correctBoundaryConditions();

	