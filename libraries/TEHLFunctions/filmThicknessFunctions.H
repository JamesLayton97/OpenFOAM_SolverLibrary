// Integral function for surface defomration performed
void pressureIntegration(volScalarField& p, volScalarField& h, volScalarField& thetaField, scalar R, scalar c, scalar E, int cellStart, int cellEnd);

// Integral function for surface defomration performed in parallel
void pressureIntegrationParallel(List<scalarField> pList, List<vectorField> pCenters, List<scalarField> pVolumeList, List<scalarField> thetaList, volScalarField& h, volScalarField& thetaField, scalar R, scalar c, scalar E);

// Intregration of pressure for single deformation point
scalar pressureIntegrate(volScalarField& p, vector C0, scalar R, scalar c, scalar E);

// Calculates the distance along the curved surface 
scalar G(scalar dx, scalar dz);
