#include "TEHLClasses.H"

convergenceTEHL::convergenceTEHL
(
 	IOdictionary dict
)
:
	relRhoTolerance_(dict.lookupOrDefault<scalar>("relRhoTolerance", 0)),
	relRhoRelTol_(dict.lookupOrDefault<scalar>("relRhoRelTol", 0)),
	pTolerance_(readScalar(dict.lookup("pTolerance"))),
	TTolerance_(readScalar(dict.lookup("TTolerance"))),
	TRegionTolerance_(readScalar(dict.lookup("TRegionTolerance"))),

	pError_(1),
	TError_(1),

	pCountMax_(readScalar(dict.lookup("pMaxIterations"))),
	TCountMax_(readScalar(dict.lookup("TMaxIterations"))),
	TEHLCountMax_(readScalar(dict.lookup("TEHLMaxIterations"))),

	pInitialResidual_(1),
	TInitialResidual_(1),	
	pInitialResidualPrev_(2),
	TInitialResidualPrev_(2),	

	TRegionInitialResidual_(1),	
	TRegionInitialResidualPrev_(2),

	loadConverged_(false),
	loadError_(1),
	loadTolerance_(0),
	loadMaxIter_(0),	
	loadTarget_(dimensionSet(1,1,-2,0,0,0,0), 0),
	currentLoad_(dimensionSet(1,1,-2,0,0,0,0), 0),

	eUpperOrig(1),
	eLowerOrig(0),
	eUpper(1),
	eLower(0),
	deltaE_(readScalar(dict.lookup("de")))
{
	deformation 		= readBool(dict.lookup("deformation"));
	temperature 		= readBool(dict.lookup("temperature"));
	viscosity 		= readBool(dict.lookup("variableViscosity"));
	load 			= readBool(dict.lookup("loadConvergence"));
	CFD 			= readBool(dict.lookup("includeCFD"));

	const word cavitationMethodInput(dict.lookup("cavitationMethod"));

	if (cavitationMethodInput == "none")
	{
		cavitationMethod = NONE;
		Info << "\tNo Cavitation method applied" << nl << endl;
	}
	else if (cavitationMethodInput == "halfSommerfeld")
	{
		cavitationMethod = HALFSOMMERFELD;
		Info << "\tHalf-Sommerfeld Cavitation Method Selected" << nl << endl;
	}
	else if (cavitationMethodInput == "elrodAdam")
	{
		cavitationMethod = ELRODADAM;
		Info << "\tElrod Adams Cavitation Method Selected" << nl << endl;
	}
	else
	{
		Info << "****\"" << cavitationMethodInput << "\"" << " is not an available cavitation method**** \n\n"
		<< "Available methods are:" 
		<< nl << "none" 
		<< nl << "halfSommerfeld"
		<< nl << "elrodAdam"
		<< endl;
	}
};

/******* Destructor ****************/

convergenceTEHL::~convergenceTEHL()
{};

/******* Get Load Values ***********/
void convergenceTEHL::setupLoadConvergence
(
	IOdictionary& dict
)
{
	loadTolerance_ 		= dict.lookupOrDefault<scalar>("tolerance", 1e-2);
	loadTarget_.value() 	= dict.lookup<scalar>("target");
	loadMaxIter_		= dict.lookupOrDefault<scalar>("maxIteration", 10);
	eUpperOrig		= dict.lookupOrDefault<scalar>("eUpper", 1);
	eLowerOrig		= dict.lookupOrDefault<scalar>("eLower", 0);

	eUpper = eUpperOrig;
	eLower = eLowerOrig;
};

void convergenceTEHL::updatePInitialResidual(scalar newResidual)
{
	pInitialResidualPrev_ = pInitialResidual_;
	pInitialResidual_= newResidual;
};

void convergenceTEHL::updateTInitialResidual(scalar newResidual)
{
	TInitialResidualPrev_ = TInitialResidual_;
	TInitialResidual_= newResidual;
};

void convergenceTEHL::updateTRegionInitialResidual(scalar newResidual)
{
	TRegionInitialResidualPrev_ = TRegionInitialResidual_;
	TRegionInitialResidual_= newResidual;
};

/******* Access Functions ***********/
scalar& convergenceTEHL::pInitialResidual()
{
	return pInitialResidual_;
};

scalar& convergenceTEHL::TInitialResidual()
{
	return TInitialResidual_;
};

scalar convergenceTEHL::TRel()
{
	scalar TRel_ = mag(TRegionInitialResidual_ - TRegionInitialResidualPrev_)/TRegionInitialResidual_;

	return TRel_;
};

scalar& convergenceTEHL::TRegionInitialResidual()
{
	return TRegionInitialResidual_;
};

scalar& convergenceTEHL::TRegionInitialResidualOld()
{
	return TRegionInitialResidualPrev_;
};

scalar& convergenceTEHL::pCountMax()
{
	return pCountMax_;
};

scalar& convergenceTEHL::TCountMax()
{
	return TCountMax_;
};

scalar& convergenceTEHL::TEHLCountMax()
{
	return TEHLCountMax_;
};

scalar& convergenceTEHL::relRhoTolerance()
{
	return relRhoTolerance_;
};

scalar& convergenceTEHL::relRhoRelTol()
{
	return relRhoRelTol_;
};

scalar& convergenceTEHL::pTolerance()
{
	return pTolerance_;
};

scalar& convergenceTEHL::TTolerance()
{
	return TTolerance_;
};

scalar& convergenceTEHL::TRegionTolerance()
{
	return TRegionTolerance_;
};

scalar& convergenceTEHL::pError()
{
	return pError_;
};

scalar& convergenceTEHL::TError()
{
	return TError_;
};

scalar& convergenceTEHL::loadTolerance()
{
	return loadTolerance_;
};

scalar& convergenceTEHL::loadMaxIter()
{
	return loadMaxIter_;
};

dimensionedScalar& convergenceTEHL::currentLoad()
{
	return currentLoad_;
};

/******* Member Functions ***********/
void convergenceTEHL::updateCurrentLoad(scalar load)
{
	currentLoad_.value() = load;
	updateLoadError();
};

void convergenceTEHL::reset()
{	
	pInitialResidual_ = 1;
	TInitialResidual_ = 1;	
	pInitialResidualPrev_ = 2;
	TInitialResidualPrev_ = 2;
	pError_ = 1;
	TError_ = 1;
};

void convergenceTEHL::resetE()
{	
	eUpper = eUpperOrig;
	eLower = eLowerOrig;
};

void convergenceTEHL::resetERel(scalar& e)
{	
	eUpper = e+deltaE_;
	eLower = e-deltaE_;
};

bool convergenceTEHL::allConverged()
{	
	bool converged = true;
	
	Info << "TInitialResidual = " << TInitialResidual_ << ", TTol = " << TTolerance_ << endl;
	Info << "pInitialResidual = " << pInitialResidual_ << ", pTol = " << pTolerance_ << endl;

	if(temperature && (TInitialResidual_ > TTolerance_))
	{
		converged = false;
	}

	if (pInitialResidual_ > pTolerance_)
	{
		converged = false;
	}

	if (converged)
	{
		Info << "TEHL Converged!!!" << endl;
	}
	
	return converged;
};

bool convergenceTEHL::TRegionConverged()
{	
	if (TRegionInitialResidual_ < TRegionTolerance_)
	{
		return true;
	}
	else
	{
		return false;
	}
};


bool convergenceTEHL::loadConverged()
{	
	if (mag(loadError_) < loadTolerance_)
	{
		return true;
	}
	
	return false;	
};

void convergenceTEHL::updateLoadError()
{
	loadError_ = ((currentLoad_ - loadTarget_)/loadTarget_).value();
	Info << "Load Error = " << loadError_ << endl;
};

scalar convergenceTEHL::loadError()
{
	return loadError_;
};

bool& convergenceTEHL::includeDeformation()
{
	return deformation;
};

bool& convergenceTEHL::includeTemperature()
{
	return temperature;
};

bool& convergenceTEHL::variableViscosity()
{
	return viscosity;
};

bool& convergenceTEHL::loadConvergence()
{
	return load;
};

bool& convergenceTEHL::includeCFD()
{
	return CFD;
};


