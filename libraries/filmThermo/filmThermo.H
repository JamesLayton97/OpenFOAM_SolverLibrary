/* Thermo calss to hold the fields in the film region for the TEHL methodology. 
 * Constructs coordinate system relative to origin and axis relating to radial, axial and circumferential directions.
 * Contains two overlapping meshes:
 * 	- 2D mesh for the pressure solution (requires 1 cell radial discretisation)
 * 	- 3D mesh for the temperature solution (requires radial discretisation)
 * Base class without Elrod Adams cavitation method
 */


//#include "rhoReactionThermo.H"
#include "fluidThermo.H"

class TEHLThermo
: public fluidThermo
{
private:
	// Film Meshes
	const fvMesh& mesh2D_;
	const fvMesh& mesh3D_;

	// Coordinate system fields
	volVectorField radialVec_;
	volVectorField cirfumferentialVec_;
	volVectorField axialVec_;
		
	// scalar Fields
	volScalarField h2D_;
	volScalarField p2D_;
	volScalarField T2D_;
	volScalarField h_;		// Film thickness
	volScalarField hGeo_;		// Film thickness from shaft eccentricity 
	volScalarField hDef_;		// Deformation field
	volScalarField theta_;	

	volVectorField U2D_;

	volScalarField g_;
	volScalarField relRho_;
	volScalarField B_;
	volScalarField mu_;
	volScalarField mu2D_;
	volScalarField rho_;
	volScalarField cp_;
	volScalarField k_;

	volVectorField Couette_;	// Components of velocity
	volVectorField Poiseuille_;
	volVectorField U_;




	/* OLD SETUP */
	// clearance of case and mesh
	dimensionedScalar c_;	
	dimensionedScalar meshC_;

	// Radius and length
	dimensionedScalar R_;	
	dimensionedScalar L_;
	dimensionedScalar omega_;	

	// Solid constants	
	dimensionedScalar E_;
	dimensionedScalar alphaShaft_;
	dimensionedScalar alphaBush_;
	dimensionedScalar TRef_;

	// Fluid constants
	dimensionedScalar BConst_;
	dimensionedScalar BCav_;
	dimensionedScalar PCav_;	
	dimensionedScalar rhoConst_;	
	dimensionedScalar rhoAirConst_;
	dimensionedScalar muConst_;
	dimensionedScalar muAirConst_;
	dimensionedScalar cpConst_;
	dimensionedScalar cpAirConst_;
	dimensionedScalar kConst_;
	dimensionedScalar kAirConst_;
	dimensionedScalar dRhoT_;
	dimensionedVector muTCoeff_;
	dimensionedVector muPCoeff_;
	
	dimensionedScalar bushT_;
	dimensionedVector misalignmentCenter_;
	dimensionedScalar hDil_;

	// non-dimensioned variables
	scalar misalignmentAngle_;	
	scalar psi_;
	scalar loadAngle_;
	scalar gFactor_;
	scalar gLimit_;
	
	int nThreads_;
	
	bool deformation;
	bool temperature;	
	bool dilation;
	bool misalignment;
	bool unwrapped;
	bool includePoiseuille;
	
	//Patch Labels
	label shaftPatchId_;
	label bushPatchId_;
	label sidesPatchId_;
	label inletPatchId_;
	label outletPatchId_;

	labelListList cellAddressing_;
	labelListListList faceAddressing_;
	wordList interpBoundaries_;

public:    	   	
	/********************* Constructor **************************/  
	TEHLThermo
	(
		const fvMesh& mesh2D,
		const fvMesh& mesh,
		IOdictionary settingDict,
		IOdictionary controlDict,
		IOdictionary solidProperties,
		IOdictionary fluidProperties
	);
	
	/********************* Destructor **************************/  
	~TEHLThermo();

	/********************* Member Functions ********************/  
	void testPatchId();
	void calculateTheta();
	void calculateRadial();
	void calculateU2D();
	void calculateU();
    	void calculateGeometric();
    	void calculateDeformation();
	void calculateDilation();
    	void calculateMu();
    	void calculateB();
	void calculateP();
    	void elrodAdamSwitch();	
	void setSwitchFunction(scalar value=-1);
	void calculateCavFields();
    	void updateFilmThickness();
    	void reset();
    	void updateVelocity();
	void updateRho();
	void updatek();
	void updatecp();
	void updateFilmFraction();
    	dimensionedScalar calculateLoadVertical();
    	dimensionedScalar calculateLoad();
	vector calculateForce();
	void updatePsi(scalar dPsi);

	void calculateAddressing();
	void interpolateT();
	void interpolateTo3D();
	void interpolatePToRelRho();

	void viscosityIntergals();

	/********************* Calculated Fields *******************/
	volScalarField hByc();

	/********************* Access Functions ********************/  	
	volScalarField& p();
	volScalarField& T();
	volScalarField& T2D();
    	volScalarField& relRho();
	volScalarField& relRho2D();
	volScalarField& p2D();
    	volScalarField& filmFraction();
    	volScalarField& g();
    	volScalarField& g3D();
	volScalarField& B();
    	volScalarField& h();
    	volScalarField& hGeo();
    	volScalarField& hDef();
	volScalarField& h2D();
    	volScalarField& theta();
    	volScalarField& rhoFilm();
    	volScalarField& muFilm();
	volScalarField& muFilm2D();
    	volScalarField& cp();
	volScalarField& F2();
	volScalarField& F3();

	surfaceScalarField massFlux();
	surfaceScalarField kSurf();

    	volVectorField& radial();
    	volVectorField& Couette();
    	volVectorField& Poiseuille();
	volVectorField& U();
	volVectorField& U2D();
	
	dimensionedScalar& PCav();
    	dimensionedScalar& c();
    	dimensionedScalar& R();
    	dimensionedScalar& L();
	dimensionedScalar& meshC();

    	dimensionedScalar& dRhoT();
	
    	dimensionedScalar& rhoConst();
    	dimensionedScalar& muConst();
    	dimensionedScalar& omega();
	dimensionedScalar UShaft();
    	dimensionedScalar& bushT();
	dimensionedScalar& hDil();
    	scalar& psi();
	scalar&loadAngle();
    	scalar e;

    	bool& includeDeformation();
	bool& includeTemperature();
	bool& includeMisalignment();
	bool& unwrappedCase();

	label& shaftPatchId();
	label& bushPatchId();
	label& sidesPatchId();
	label& inletPatchId();
	label& outletPatchId();
	
	
	/********************* Derived Functions ********************/
	virtual void correct();
	
	virtual word thermoName() const;
	
	virtual bool incompressible() const;
	
	virtual bool isochoric() const;
	
	virtual tmp<volScalarField> rho() const;
	
	virtual tmp<scalarField> rho(const label patchi) const;
	
	virtual volScalarField& he();
	
	virtual tmp<volScalarField> he(const volScalarField& p, const volScalarField& T) const;

	virtual tmp<scalarField> he(const scalarField& T,const labelList& cells) const;

	virtual tmp<scalarField> he(const scalarField& T,const label patchi) const;

	virtual const volScalarField& he() const;
	
	virtual tmp<volScalarField> hs() const;

	virtual tmp<volScalarField> hs(const volScalarField& p,const volScalarField& T) const;

	virtual tmp<scalarField> hs(const scalarField& T,const labelList& cells) const;

	virtual tmp<scalarField> hs(const scalarField& T,const label patchi) const;
	
	virtual tmp<volScalarField> ha() const;

	virtual tmp<volScalarField> ha(const volScalarField& p,const volScalarField& T) const;

	virtual tmp<scalarField> ha(const scalarField& T,const labelList& cells) const;

	virtual tmp<scalarField> ha(const scalarField& T,const label patchi) const;
	
	virtual tmp<scalarField> THE(const scalarField& h,const scalarField& T0,const labelList& cells) const;

	virtual tmp<scalarField> THE(const scalarField& h, const scalarField& T0,const label patchi) const;
	
	virtual tmp<volScalarField> hc() const;
	
	virtual tmp<volScalarField> Cp() const;
	
	virtual tmp<scalarField> Cp(const scalarField& T,const label patchi) const;

	virtual tmp<volScalarField> Cv() const;

	virtual tmp<scalarField> Cv(const scalarField& T,const label patchi) const;

	virtual tmp<volScalarField> Cpv() const;

	virtual tmp<scalarField> Cpv(const scalarField& T,const label patchi) const;

	virtual tmp<volScalarField> CpByCpv() const;

	virtual tmp<scalarField> CpByCpv(const scalarField& T,const label patchi) const;

	virtual tmp<volScalarField> kappa() const;

	virtual tmp<scalarField> kappa(const label patchi) const;
	
	virtual tmp<volScalarField> alphahe() const;
	
	virtual tmp<scalarField> alphahe(const label patchi) const;
	
	virtual void correctRho(const volScalarField& deltaRho);
	
	virtual const volScalarField& psi() const;
	
	virtual tmp<volScalarField> gamma() const;

	virtual tmp<scalarField> gamma(const scalarField& T,const label patchi) const;
	
	virtual tmp<volScalarField> W() const;
	
	virtual tmp<scalarField> W(const label patchi) const;
	
	virtual tmp<volScalarField> mu() const;
	
	virtual tmp<scalarField> mu(const label patchi) const;
	
	virtual tmp<volScalarField> kappaEff(const volScalarField&) const;

	virtual tmp<scalarField> kappaEff(const scalarField& alphat,const label patchi) const;

	virtual tmp<volScalarField> alphaEff(const volScalarField& alphat) const;

	virtual tmp<scalarField> alphaEff(const scalarField& alphat,const label patchi) const;

};
