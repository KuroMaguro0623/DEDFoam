
m_vap = psat*Foam::sqrt(1e-3*mw/(2*Foam::constant::mathematical::pi*R*T));

Q_vapor = -m_vap*L_vapor*mag(fvc::grad(alpha1))*factor;
Q_vapor.correctBoundaryConditions();

// Info << "min/max(p_sat) = " << gMin(psat) << ", " << gMax(psat) << endl;
// Info << "min/max(m_vap) = " << gMin(m_vap) << ", " << gMax(m_vap) << endl;
// Info << "min/max(Q_vapor) = " << gMin(Q_vapor) << ", " << gMax(Q_vapor) << endl;
Info << "++Vapour heat has been calculated." << endl;
