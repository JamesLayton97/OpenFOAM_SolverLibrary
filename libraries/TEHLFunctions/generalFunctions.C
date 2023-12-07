#include "TEHLFunctions.H"

scalar convergeTest(volScalarField& previous, volScalarField &current, bool& converged, scalar& condition)
{
	scalar error(0);

	forAll(current, i)
	{     
		error += Foam::sqrt(Foam::pow(current[i] - previous[i], 2));
	}

	error = error/current.mesh().nCells()/max(current).value();

	if (error < condition)
	{
		converged = true;
	}
	else
	{
		converged = false;
	}

	return error;
}

scalar calculateFieldError(const volScalarField& previous, volScalarField& current)
{
	scalar error(0);

	forAll(current, i)
	{     
		error += Foam::mag(current[i] - previous[i])/current[i];
	}

	scalar nCells = current.mesh().nCells();
	scalar fieldMax = max(current).value();

	if (Pstream::parRun())
	{
		reduce(error, sumOp<scalar>());
		reduce(nCells, sumOp<scalar>());
		reduce(fieldMax, maxOp<scalar>());
	}

	error = error/(nCells);

	return error;
}

dimensionedScalar calculateTorque(volScalarField& mu, volScalarField& h, volVectorField& U, dimensionedScalar& c, dimensionedScalar& R)
{
	scalar Area;
	dimensionedScalar torque(dimensionSet(1,2,-2,0,0,0,0), 0);

	forAll(mu, cellID)
	{
		Area 		= mu.mesh().V()[cellID]/c.value(); 
		torque.value() 	+= mu[cellID]*Area*R.value()*(mag(U[cellID])/h[cellID]);
	}

	return torque;
}

void interpolateToBoundaries(volScalarField& field)
{
	forAll(field.mesh().boundary(), patchId)
	{
		const scalarField fieldInternal = field.boundaryField()[patchId].patchInternalField();

		forAll(field.mesh().boundary()[patchId], faceId)
		{
			field.boundaryFieldRef()[patchId][faceId] = fieldInternal[faceId];
		}
	}
}

label nearestBoundaryFace(const vector& faceCenter, const vectorField& boundary)
{
	label closestFaceId(0);
	scalar closestDistance = mag(faceCenter - boundary[0]);

	forAll(boundary, faceId)
	{
		scalar faceDistance = mag(faceCenter - boundary[faceId]);

		if (faceDistance < closestDistance)
		{
			closestFaceId = faceId;
			closestDistance = faceDistance;
		}
	}
	
	return closestFaceId;	
}
