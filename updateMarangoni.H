const dimensionedScalar deltaN("deltaN", dimless/dimLength, 1e-6);

volVectorField gradT(fvc::grad(T));
// surfaceVectorField gradTSf(fvc::interpolate(gradT));

// volVectorField gradAlpha(fvc::grad(alpha1));
// surfaceVectorField gradAlphaSf(fvc::interpolate(gradAlpha));

// surfaceVectorField nHatSf(gradAlphaSf/(mag(gradAlphaSf) + deltaN));
// surfaceScalarField gradTdotn(nHatSf & gradTSf);

// surfaceScalarField MarangoniSf = dSigma_dT*(fvc::snGrad(T)*mag(gradAlphaSf) - fvc::snGrad(alpha1)*gradTdotn);
// Marangoni = fvc::reconstruct(MarangoniSf*mesh.magSf());


volVectorField gradAlpha = fvc::grad(alpha1);
volVectorField nHat(gradAlpha/(mag(gradAlpha) + deltaN));

Marangoni = dSigma_dT*(gradT - nHat*(nHat&gradT))*mag(gradAlpha)*factor;
Marangoni.correctBoundaryConditions();
