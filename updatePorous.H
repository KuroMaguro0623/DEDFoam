const volVectorField& C = mesh.C();


forAll(C, i)
{
	if( T[i]<=Ts.value())
	{
		beta[i] = 0 ;
	}
	else if( T[i]>=Tl.value())
	{
		beta[i] = 1 ;
	}
	else 
		beta[i] = ( T[i]-Ts.value() )/( Tl.value()-Ts.value() ) ;
	
	if( alpha1[i] <= 0.01)
		beta[i] = 1;
}
beta.correctBoundaryConditions();
A_c = -A_mush*Foam::sqr(1 - beta)/(Foam::pow(beta, 3) + bb); //implicit's written

