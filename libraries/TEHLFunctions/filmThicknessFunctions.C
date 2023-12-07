#include "TEHLFunctions.H"

// Deformation integral calculation in the case of a serial run with multithreading
void pressureIntegration(volScalarField& p, volScalarField& h, volScalarField& thetaField, scalar R, scalar c, scalar E, int cellStart, int cellEnd)
{
	// Loop through the cell range defined in the parameters
	for(int cellId = cellStart; cellId < cellEnd; cellId++)
	{
		// Initialise variable for calculating the sum of the influence from each cell to the deformation of the current cell
	    	scalar sum(0);
		
		// Store addresses of the current cell position vector and angular position
		const vector& cell0Position = h.mesh().C()[cellId];
		scalar& theta0 = thetaField[cellId];
		
		// Loop through all cells for contribution to deformation of cell of interest
		forAll(p, i)
		{    
			// define values for contributing cell
			const vector& cellPosition = p.mesh().C()[i];
			scalar& theta 	= thetaField[i];

			// Calculate and angular distance between the cells
			scalar dTheta 	= Foam::sqrt(Foam::pow(theta0 - theta, 2));

			// If statement to get the smallest difference
			if(dTheta > 3.142)
			{
				dTheta = 6.284 - dTheta;
			}

			// Only include cells within 1 quarter rotation, assumes force vectors angled greater >90 degrees from the face normal do not contribute
			//if(dTheta < 1.57)
			//{
				// Calculate the distat
				scalar dx = R*dTheta;
				scalar dz = cell0Position.z() - cellPosition.z();
				scalar dR = Foam::sqrt(Foam::pow(dx, 2) + Foam::pow(dz, 2));	

				// if dx and dz are 0, the denominator becomes 0 in the intergral, therefore a minimum value is applied to prevent a floating point error
				if (dR < 1e-5)
				{
					dR = 1e-5;
				}

				// calculate the planar surface area of the cell 
				scalar A = p.mesh().V()[i]/c;

				// apply gauss centered rectangular method of integration
				sum += (p[i] * A)/dR;
			//}
		}  

		// Calculate the deformation of the cell
		h[cellId] = (2/(3.14159*E)) * sum;
	}
}

// Applies the deformation integral function for a case in parallel
void pressureIntegrationParallel(List<scalarField> pList, List<vectorField> pCentersList, List<scalarField> pVolumeList, List<scalarField> thetaList, volScalarField& h, volScalarField& thetaField, scalar R, scalar c, scalar E)
{
	// Loop through all cells in deformation field
	forAll(h, cellId)
	{
		// Initialise variable to store the sum of the influence of all the cells 
	    	scalar sum(0);
		
		// Get the address of the current cells position vector and angular position
		const vector& cell0Position = h.mesh().C()[cellId];
		scalar& theta0 = thetaField[cellId];
		
		// Loop throught all the processors to get the contribution from all cells across the domain
		forAll(pList, proc)
		{    
			// Get the lists of information for the cells on a processor
			scalarField& p = pList[proc];
			vectorField& pCenters = pCentersList[proc];
			scalarField& pVolumes = pVolumeList[proc];
			scalarField& thetaVals = thetaList[proc];

			// Loop all cells on the processor to calculate the influence of each cell 
			forAll(p, pId)
			{
				// Get the address of the cell position vector
				vector& cellPosition = pCenters[pId];

				// Estimate the area of the cell against the shaft/bush surfaces (average between the two)
				scalar A = pVolumes[pId]/c;

				// Get the address of the cells angular position
				scalar& theta = thetaVals[pId];

				// Calculate the shortest angular distance between the cells
				scalar dTheta 	= Foam::sqrt(Foam::pow(theta0 - theta, 2));
				if (dTheta > 3.142)
				{
					dTheta = 6.284 - dTheta;
				}

				// Calculate the distance along the curved surface between the cells
				scalar dx = R*dTheta;
				scalar dz = cell0Position.component(2) - cellPosition.component(2);
				scalar dR = Foam::sqrt(Foam::pow(dx, 2) + Foam::pow(dz, 2));

				// if dx and dz are 0, the denominator becomes 0 in the intergral, therefore a minimum value is applied to prevent a floating point error
				if (dR < 1e-5)
				{
					dR = 1e-5;
				}
				
				// Apply gauss centered rectangular method of integration
				sum += (p[pId] * A)/dR;
			}
		}  
		
		// Calcualte the deformation of the cell
		h[cellId] = (2/(3.14159*E)) * sum;
	}
}

scalar pressureIntegrate(volScalarField& p, vector C0, scalar R, scalar c, scalar E)
{
	scalar pSum(0);
	scalar theta0 = Foam::atan2(C0.x(), -1*C0.y()); 
	
	// Loop through all cells for contribution to deformation of cell of interest
	forAll(p, pCellId)
	{    
		// define values for contributing cell
		const vector& C = p.mesh().C()[pCellId];
			
		scalar theta 	= Foam::atan2(C.x(), -1*C.y());
		scalar dTheta 	= Foam::sqrt(Foam::pow(theta0 - theta, 2));

		if (dTheta > 3.142)
		{
			dTheta = 6.284 - dTheta;
		}

		scalar dx = R*dTheta;
		scalar dz = C0.z() - C.z();
		

		// calculate the average radial normal surface area of the cell 
		scalar A = p.mesh().V()[pCellId]/c;

		// apply gauss centered rectangular method of integration
		pSum += p[pCellId] * A * G(dx, dz);
	}  

	pSum *= (2/(3.14159*E));

	return pSum;
}

scalar G(scalar dx, scalar dz)
{
	scalar coeff, distance;

	distance = Foam::sqrt(Foam::pow(dx, 2) + Foam::pow(dz, 2));

	if (distance > 1e-6)
	{
		coeff = 1/distance;
	}
	else
	{
		coeff = 0;
	}

	return coeff;
}
