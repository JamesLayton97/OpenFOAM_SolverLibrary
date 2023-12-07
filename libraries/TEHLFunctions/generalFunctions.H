scalar convergeTest(volScalarField& previous, volScalarField &current, bool& converged, scalar& condition);

scalar calculateFieldError(const volScalarField& previous, volScalarField &current);

dimensionedScalar calculateLoadRig(volScalarField& p, volScalarField& Theta, dimensionedScalar& c);

dimensionedScalar calculateLoadJournal(volScalarField& p, volScalarField& Theta, dimensionedScalar& c);

dimensionedScalar calculateTorque(volScalarField& mu, volScalarField& h, volVectorField& U, dimensionedScalar& c, dimensionedScalar& R);

void interpolateToBoundaries(volScalarField& field);

label nearestBoundaryFace(const vector& faceCenter, const vectorField& boundary);
