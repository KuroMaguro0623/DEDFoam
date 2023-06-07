#include "updatePorous.H"
#include "updateMarangoni.H"
#include "updateRecoil.H"

volVectorField surfaceTension = fluid.surfaceTensionCoeff(phasePairKey(fluid.phase(0).name(), fluid.phase(1).name()))*
								fluid.K(alpha1, alpha2)*
								fvc::grad(alpha1)*factor;

fvVectorMatrix UEqn
(
    fvm::ddt(rho, U)
  + fvm::div(rhoPhi, U)
//- fvm::Sp(fvc::ddt(rho) + fvc::div(rhoPhi), U)
  + turbulence->divDevRhoReff(U)
  ==
    fvOptions(rho, U)
  + fvm::Sp(A_c, U)
  + 0.54*recoil
  + Marangoni
  + surfaceTension
);

UEqn.relax();

fluid.addInterfacePorosity(UEqn);

if (pimple.momentumPredictor())
{
    solve
    (
        UEqn
      ==
        fvc::reconstruct
        (
            (
                // fluid.surfaceTensionForce()
				// + MarangoniSf
                - ghf*fvc::snGrad(rho)
                - fvc::snGrad(p_rgh)
            ) * mesh.magSf()
        )
    );

    fvOptions.correct(U);
    K = 0.5*magSqr(U);
}
