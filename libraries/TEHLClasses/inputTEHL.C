#include "TEHLClasses.H"

inputTEHL::inputTEHL
(
	fvMesh& filmMesh,
	Time& runTime
):
	settingsDict
	(
	    IOobject
	    (
		"settings",
		runTime.constant(),
		filmMesh,
		IOobject::MUST_READ,
		IOobject::NO_WRITE
	    )
	), 
	
	filmControlDict
	(
	    IOobject
	    (
		"solverControl",
		runTime.system(),
		filmMesh,
		IOobject::MUST_READ,
		IOobject::NO_WRITE
	    )
	),

	solidConstantsDict
	(
	    IOobject
	    (
		"solidConstants",
		runTime.constant(),
		filmMesh,
		IOobject::MUST_READ,
		IOobject::NO_WRITE
	    )
	), 

	fluidConstantsDict
	(
	    IOobject
	    (
		"fluidConstants",
		runTime.constant(),
		filmMesh,
		IOobject::MUST_READ,
		IOobject::NO_WRITE
	    )
	),

	outputControlDict
	(
	    IOobject
	    (
		"outputControl",
		runTime.system(),
		filmMesh,
		IOobject::MUST_READ,
		IOobject::NO_WRITE
	    )
	)
{
};

inputTEHL::~inputTEHL()
{};

IOdictionary inputTEHL::settings()
{
	return settingsDict;
};

IOdictionary inputTEHL::filmControl()
{
	return filmControlDict;
};

IOdictionary inputTEHL::solidConstants()
{
	return solidConstantsDict;
};

IOdictionary inputTEHL::fluidConstants()
{
	return fluidConstantsDict;
};

IOdictionary inputTEHL::outputControl()
{
	return outputControlDict;
};



