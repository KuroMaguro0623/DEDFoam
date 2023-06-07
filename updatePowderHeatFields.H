const volScalarField& dmdt = mesh.lookupObject<volScalarField>("powderDistribution");

Q_powder =  dmdt* (h_powder - rhoCp/rho*T)*pos(T - Tl);
Q_powder.correctBoundaryConditions();
 
Info << "min/max(Q_powder) = " << gMin(Q_powder) << ", " << gMax(Q_powder) << endl;
Info << "++Powder heat has been calculated." << endl;