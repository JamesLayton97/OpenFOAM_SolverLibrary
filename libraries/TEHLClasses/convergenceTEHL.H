class convergenceTEHL
{
	private:

	enum cavitationMethodType 
	{
		NONE,
		HALFSOMMERFELD,
		ELRODADAM,
	};
	bool deformation;
	bool temperature;
	bool viscosity;
	bool load;
	bool CFD;
	cavitationMethodType cavitationMethod;

	// Input tolerance values to individual fields
	scalar relRhoTolerance_;
	scalar relRhoRelTol_;
	scalar pTolerance_;
	scalar TTolerance_;
	scalar TRegionTolerance_;

	scalar pError_;
	scalar TError_;

	scalar pCountMax_;
	scalar TCountMax_;
	scalar TEHLCountMax_;

	scalar pInitialResidual_;
	scalar TInitialResidual_;	
	scalar pInitialResidualPrev_;
	scalar TInitialResidualPrev_;
	scalar TRegionInitialResidual_;	
	scalar TRegionInitialResidualPrev_;

	// Boolean variables to indicate converged solutions
	bool TRegionConverged_;	
	
	// Variables for if program is seeking input load
	bool loadConverged_;
	scalar loadError_;
	scalar loadTolerance_;
	scalar loadMaxIter_;
	dimensionedScalar loadTarget_;
	dimensionedScalar currentLoad_;

	public:

	scalar eUpperOrig;
	scalar eLowerOrig;

	scalar eUpper;
	scalar eLower;
	scalar deltaE_;

	/******* Constructor *************/
	convergenceTEHL
	(
	 	IOdictionary dict
	);
	
	/******* Destructor *************/

	~convergenceTEHL();

	/******* Get Load Values ***********/
	void setupLoadConvergence(IOdictionary& dict);

	/******* Update Residual Values ******/
	void updatePInitialResidual(scalar newResidual);
	void updateTInitialResidual(scalar newResidual);
	void updateTRegionInitialResidual(scalar newResidual);

	/******* Access Functions ***********/
	scalar& pInitialResidual();
	scalar& TInitialResidual();
	scalar& TRegionInitialResidual();
	scalar& TRegionInitialResidualOld();
	scalar TRel();
	scalar& pCountMax();
	scalar& TCountMax();
	scalar& TEHLCountMax();
	scalar& relRhoTolerance();
	scalar& relRhoRelTol();
	scalar& pTolerance();
	scalar& TTolerance();
	scalar& TRegionTolerance();
	scalar& pError();
	scalar& TError();
	scalar& loadTolerance();
	scalar& loadMaxIter();

	dimensionedScalar& currentLoad();

	bool& includeDeformation();
	bool& includeTemperature();
	bool& variableViscosity();
	bool& loadConvergence();
	bool& includeMisalignment();
	bool& includeCFD();
	bool& varyShaftTemp();
	cavitationMethodType getCavitationMethod()
	{
		return cavitationMethod;
	};

	/******* Member Functions ***********/

	void updateCurrentLoad(scalar load);
	void reset();
	void resetE();
	void resetERel(scalar& e);
	
	bool allConverged();
	bool TRegionConverged();
	
	bool loadConverged();
	void updateLoadError();
	scalar loadError();
};
