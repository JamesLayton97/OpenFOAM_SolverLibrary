#include "TEHLClasses.H"
#include "meshSearch.H"
#include <thread>

/*************************************Constructor****************************************/
TEHLThermo::TEHLThermo
(
	const fvMesh& mesh2D,
	const fvMesh& mesh3D,
	IOdictionary settingDict,
	IOdictionary controlDict,
	IOdictionary solidProperties,
	IOdictionary fluidProperties
)
:
	fluidThermo(mesh3D, ""),

	c_("c", dimensionSet(0,1,0,0,0,0,0), settingDict),
	meshC_("meshC", dimensionSet(0,1,0,0,0,0,0), settingDict),
	R_("R", dimensionSet(0,1,0,0,0,0,0), settingDict),
	L_("L", dimensionSet(0,1,0,0,0,0,0), settingDict),
	omega_("omega", dimensionSet(0,0,-1,0,0,0,0), settingDict),
	E_(dimensionSet(1,-1,-2,0,0,0,0), 0),
	alphaShaft_("alphaShaft", dimensionSet(0,0,0,-1,0,0,0), solidProperties),
	alphaBush_("alphaBush", dimensionSet(0,0,0,-1,0,0,0), solidProperties),
	TRef_("Tref", dimensionSet(0,0,0,1,0,0,0), solidProperties),
	BConst_("BConst", dimensionSet(1,-1,-2,0,0,0,0), fluidProperties),
	BCav_("BCav", dimensionSet(1,-1,-2,0,0,0,0), fluidProperties),
	PCav_("PCav", dimensionSet(1,-1,-2,0,0,0,0), fluidProperties),
	rhoConst_("rhoConst", dimensionSet(1,-3,0,0,0,0,0), fluidProperties),
	rhoAirConst_("rhoAirConst", dimensionSet(1,-3,0,0,0,0,0), fluidProperties),
	muConst_("muConst", dimensionSet(1,-1,-1,0,0,0,0), fluidProperties),
	muAirConst_("muAirConst", dimensionSet(1,-1,-1,0,0,0,0), fluidProperties),
	cpConst_("cpConst", dimensionSet(0,2,-2,-1,0,0,0), fluidProperties),
	cpAirConst_("cpAirConst", dimensionSet(0,2,-2,-1,0,0,0), fluidProperties),
	kConst_("kConst", dimensionSet(1,1,-3,-1,0,0,0), fluidProperties),
	kAirConst_("kAirConst", dimensionSet(1,1,-3,-1,0,0,0), fluidProperties),	
	dRhoT_("dRhoT", dimensionSet(1,-3,0,-1,0,0,0), fluidProperties),	
	muTCoeff_("muTCoeff", dimless, fluidProperties),
	muPCoeff_("muPCoeff", dimless, fluidProperties),
	bushT_(dimensionSet(0,0,0,1,0,0,0), 300),
	misalignmentCenter_(dimensionSet(0,1,0,0,0,0,0), vector(0,0,0)),	
	hDil_(dimensionSet(0,1,0,0,0,0,0), 0),	

	mesh2D_(mesh2D),
	mesh3D_(mesh3D),

	T2D_
	(
	    IOobject
	    (
		"T",
		mesh2D.time().timeName(),
		mesh2D,
		IOobject::MUST_READ,
		IOobject::AUTO_WRITE
	    ),
	  mesh2D 
	),

	g_
	(
	    IOobject
	    (
		"g",
		mesh2D.time().timeName(),
		mesh2D,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	    ),
	    mesh2D,
	    dimensionedScalar(dimless, 1)
	),

	g3D_
	(
	    IOobject
	    (
		"g",
		mesh3D.time().timeName(),
		mesh3D,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	    ),
	    mesh3D,
	    dimensionedScalar(dimless, 1)
	),

	relRho_
	(
	    IOobject
	    (
		"relRho",
		mesh3D.time().timeName(),
		mesh3D,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	    ),
	    mesh3D,
	    dimensionedScalar(dimless, 1)
	),

	relRho2D_
	(
	    IOobject
	    (
		"relRho",
		mesh2D.time().timeName(),
		mesh2D,
		IOobject::MUST_READ,
		IOobject::NO_WRITE
	    ),
	  mesh2D 
	),

	p2D_
	(
	    IOobject
	    (
		"p",
		mesh2D.time().timeName(),
		mesh2D,
		IOobject::MUST_READ,
		IOobject::AUTO_WRITE
	    ),
	  mesh2D 
	    //dimensionedScalar(dimensionSet(1,-1,-2,0,0,0,0), 0)
	),

	filmFraction_
	(
	    IOobject
	    (
		"filmFraction",
		mesh3D.time().timeName(),
		mesh3D,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	    ),
	    mesh3D,
	    dimensionedScalar(dimless, 1)
	),
	
	B_
	(
	    IOobject
	    (
		"B",
		mesh2D.time().timeName(),
		mesh2D,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	    ),
	    mesh2D,
	    dimensionedScalar("BConst", dimensionSet(1,-1,-2,0,0,0,0), fluidProperties)
	),

	h_
	(
	    IOobject
	    (
		"filmThickness",
		mesh3D.time().timeName(),
		mesh3D,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	    ),
	    mesh3D,
	    meshC_
	),
		
	hGeo_
	(
	    IOobject
	    (
		"rigidFilmThickness",
		mesh2D.time().timeName(),
		mesh2D,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	    ),
	    mesh2D,
	    dimensionedScalar(dimensionSet(0,1,0,0,0,0,0), 0)
	),
	
	hDef_
	(
	    IOobject
	    (
		"deformation",
		mesh2D.time().timeName(),
		mesh2D,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	    ),
	    mesh2D,
	    dimensionedScalar(dimensionSet(0,1,0,0,0,0,0), 0)
	),
	
	h2D_
	(
	    IOobject
	    (
		"filmThickness",
		mesh2D.time().timeName(),
		mesh2D,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	    ),
	    mesh2D,
	    dimensionedScalar(dimensionSet(0,1,0,0,0,0,0), 0)
	),

	theta_
	(
	    IOobject
	    (
		"theta",
		mesh2D.time().timeName(),
		mesh2D,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	    ),
	    mesh2D,
	    dimensionedScalar(dimless, 0)
	),
	
	mu_
	(
	    IOobject
	    (
		"mu",
		mesh3D.time().timeName(),
		mesh3D,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	    ),
	    mesh3D,
	    muConst_
	),
	
	mu2D_
	(
	    IOobject
	    (
		"mu",
		mesh2D.time().timeName(),
		mesh2D,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	    ),
	    mesh2D,
	    muConst_
	),
	
	rho_
	(
	    IOobject
	    (
		"rho",
		mesh3D.time().timeName(),
		mesh3D,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	    ),
	    mesh3D,
	    rhoConst_
	),

	cp_	
	(
	    IOobject
	    (
		"cp",
		mesh3D.time().timeName(),
		mesh3D,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	    ),
	    mesh3D,
	    cpConst_
	),
	
	k_
	(
	    IOobject
	    (
		"k",
		mesh3D.time().timeName(),
		mesh3D,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	    ),
	    mesh3D,
	    kConst_
	),

	F2_
	(
	    IOobject
	    (
		"F2",
		mesh2D.time().timeName(),
		mesh2D,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	    ),
	    mesh2D,
	    dimensionedScalar(dimensionSet(-1,4,1,0,0,0,0), 0)
	),

	F3_
	(
	    IOobject
	    (
		"F3",
		mesh2D.time().timeName(),
		mesh2D,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	    ),
	    mesh2D,
	    dimensionedScalar(dimLength, 0)
	),

	radial_
	(
	    IOobject
	    (
		"radial",
		mesh3D.time().timeName(),
		mesh3D,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	    ),
	    mesh3D,
	    dimensionedVector(dimless, vector(0,1,0))
	),

	Couette_
	(
	    IOobject
	    (
		"Couette",
		mesh3D.time().timeName(),
		mesh3D,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	    ),
	    mesh3D,
	    dimensionedVector(dimensionSet(0,1,-1,0,0,0,0), vector(0,0,0))
	),

	Poiseuille_
	(
	    IOobject
	    (
		"Poiseuille",
		mesh3D.time().timeName(),
		mesh3D,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	    ),
	    mesh3D,
	    dimensionedVector(dimensionSet(0,1,-1,0,0,0,0), vector(0,0,0))
	),
	
	U_
	(
	    IOobject
	    (
		"U",
		mesh3D.time().timeName(),
		mesh3D,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	    ),
	    mesh3D,
	    dimensionedVector(dimensionSet(0,1,-1,0,0,0,0), vector(0,0,0))
	),

	U2D_
	(
	    IOobject
	    (
		"U",
		mesh2D.time().timeName(),
		mesh2D,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	    ),
	    mesh2D,
	    dimensionedVector(dimensionSet(0,1,-1,0,0,0,0), vector(0,0,0))
	),

	psi_(settingDict.lookupOrDefault<scalar>("psi", 0)),
	gFactor_(settingDict.lookupOrDefault<scalar>("gFactor", 0.8)),
	gLimit_(settingDict.lookupOrDefault<scalar>("gLimit", 1e-5)),
	nThreads_(controlDict.lookupOrDefault<scalar>("nDeformationThreads", 1)),
	
	deformation(controlDict.lookupOrDefault<bool>("deformation", false)),
	temperature(controlDict.lookupOrDefault<bool>("temperature", false)),
	dilation(controlDict.lookupOrDefault<bool>("thermalDilation", false)),
	misalignment(controlDict.lookupOrDefault<bool>("misalignment", false)),
	unwrapped(controlDict.lookupOrDefault<bool>("unwrapped", false)),
	includePoiseuille(controlDict.lookupOrDefault<bool>("includePoiseuille", true)),
	
	shaftPatchId_(mesh3D.boundaryMesh().findPatchID("film_to_shaft")),
	bushPatchId_(mesh3D.boundaryMesh().findPatchID("film_to_bush")),
	sidesPatchId_(mesh3D.boundaryMesh().findPatchID("filmsides")),
	inletPatchId_(mesh3D.boundaryMesh().findPatchID("filminlet")),
	outletPatchId_(mesh3D.boundaryMesh().findPatchID("filmoutlet")),

	cellAddressing_(mesh2D.nCells()),
	faceAddressing_(mesh2D.boundaryMesh().size()),
	interpBoundaries_(controlDict.lookupOrDefault<wordList>("boundaries", wordList(0))),

	e(readScalar(settingDict.lookup("e")))
{
	Info << "Testing Patch IDs..." << endl;
	testPatchId();
	Info << "Calculating Theta..." << endl;
	calculateTheta();
	Info << "Calculating Film Heights..." << endl;
	calculateGeometric();
	Info << "Calculating Average Couette..." << endl;
	calculateU2D();
	Info << "Calculating Full Couette..." << endl;
	calculateU();
	Info << "Calculating Addressing Between 2D And 3D Meshes..." << endl;
	calculateAddressing();
	Info << "Calculating Radial Vector Field..." << endl;
	calculateRadial();

	if(misalignment)
	{
		Info << "Getting misalignment variables" << endl;
		misalignmentAngle_ = readScalar(settingDict.lookup("misalignmentAngle"));
		misalignmentCenter_ = settingDict.lookup("misalignmentCenter");
	}
	else
	{
		misalignmentAngle_ = 0;
	}
	
	scalar EShaft	= readScalar(solidProperties.lookup("EShaft"));
	scalar EBush 	= readScalar(solidProperties.lookup("EBush"));
	scalar vShaft	= readScalar(solidProperties.lookup("vShaft"));
	scalar vBush 	= readScalar(solidProperties.lookup("vBush"));

	E_.value() = 2/((1-Foam::pow(vShaft, 2))/EShaft + (1-Foam::pow(vBush, 2))/EBush);
	Info << "Effective Youngs Modulus: " << E_.value() << endl;
	
	h2D_ = hGeo_;
};

/*************************************Destructor****************************************/

TEHLThermo::~TEHLThermo()
{};

/*************************************Calculation Functions****************************************/

// Test if patch name are in the case
void TEHLThermo::testPatchId()
{
	if(shaftPatchId_ == -1)
	{
		Info << "\tNo patch named \"film_to_shaft\"" << endl;
	}
	if(bushPatchId_ == -1)
	{
		Info << "\tNo patch named \"film_to_bush\"" << endl;
	}
	if(sidesPatchId_ == -1)
	{
		Info << "\tNo patch named \"filmsides\"" << endl;
	}
	if(inletPatchId_ == -1)
	{
		Info << "\tNo patch named \"filminlet\"" << endl;
	}
	if(outletPatchId_ == -1)
	{
		Info << "\tNo patch named \"filmoutlet\"" << endl;
	}
};

// Calculate theta values in a scalar field
void TEHLThermo::calculateTheta()
{
	forAll(theta_, cellId)
	{
		if (unwrapped)
		{
			/* Unwrapped Journal Bearin */
			theta_[cellId] = theta_.mesh().C()[cellId].x()/R_.value();
		}
		else
		{
			/* Journal Bearing Geometry */
			theta_[cellId] = Foam::atan2(theta_.mesh().C()[cellId].x(), -1*theta_.mesh().C()[cellId].y());            
		}
	}

	forAll(theta_.mesh().boundary(), patchID)
	{
		forAll(theta_.mesh().boundary()[patchID], i)
		{
			if (unwrapped)
			{
				/* Unwrapped Journal Bearing */	
				theta_.boundaryFieldRef()[patchID][i] = theta_.mesh().C().boundaryField()[patchID][i].x()/R_.value();  
			}
			else
			{
				/* Journal Bearing Geometry */
				theta_.boundaryFieldRef()[patchID][i] = Foam::atan2(theta_.mesh().C().boundaryField()[patchID][i].x(), -1*theta_.mesh().C().boundaryField()[patchID][i].y());        
			}
		}
	}
};

// Calculate Radial Vecotr Field for Thermal Conductivity adjustments
void TEHLThermo::calculateRadial()
{
	if (!unwrapped)
	{
		forAll(radial_, cellId)
		{
			scalar cellx = mesh3D_.C()[cellId].x();
			scalar celly = mesh3D_.C()[cellId].y();
			scalar R = Foam::sqrt(cellx*cellx + celly*celly);

			radial_[cellId] = vector(-1*cellx/R, -1*celly/R, 0);
		}

		forAll(radial_.mesh().boundary(), patchId)
		{
			forAll(radial_.mesh().boundary()[patchId], faceId)
			{
				scalar facex = mesh3D_.Cf().boundaryField()[patchId][faceId].x();
				scalar facey = mesh3D_.Cf().boundaryField()[patchId][faceId].y();
				scalar R = Foam::sqrt(facex*facex + facey*facey);

				radial_.boundaryFieldRef()[patchId][faceId] = vector(-1*facex/R, -1*facey/R, 0);
			}
		}
	}
};


// Calculate initial velocity field
void TEHLThermo::calculateU2D()
{
	scalar x, y;
    	scalar UMag = (omega_*R_).value()/2;
    	
	forAll(theta_, cellId)
	{  
		x = Foam::cos(theta_[cellId]);
		y = Foam::sin(theta_[cellId]);

		if(unwrapped)
		{
			/* Unwrapped Journal Bearing */
			U2D_[cellId] = UMag*vector(1, 0, 0);
		}
		else
		{
			/* Journal Bearing Geometry */
			U2D_[cellId] = UMag*vector(x, y, 0);
		}
	}

	forAll(theta_.mesh().boundary(), patchID)
	{
		forAll(theta_.mesh().boundary()[patchID], cellId)
		{       
			x = Foam::cos(theta_.boundaryFieldRef()[patchID][cellId]);
			y = Foam::sin(theta_.boundaryFieldRef()[patchID][cellId]);

			if(unwrapped)
			{
				/* Unwrapped Journal Bearing */
				U2D_.boundaryFieldRef()[patchID][cellId] = UMag*vector(1, 0, 0);
			}
			else
			{
				/* Journal Bearing Geometry */
				U2D_.boundaryFieldRef()[patchID][cellId] = UMag*vector(x, y, 0);
			}
		}
	}
};

void TEHLThermo::calculateU()
{
	scalar UMag = (omega_*R_).value();
	scalar R = R_.value();
	scalar m = -1*UMag/meshC_.value();

	forAll(Couette_, cellId)
	{  
		if(unwrapped)
		{
			/* Unwrapped Journal Bearing */
			scalar y = Couette_.mesh().C()[cellId].y();
			Couette_[cellId] = ((UMag/2)-m*y)*vector(1, 0, 0);
		}
		else
		{
			/* Journal Bearing Geometry */
			scalar y = Couette_.mesh().C()[cellId].y();
			scalar x = Couette_.mesh().C()[cellId].x();
			scalar cellTheta = Foam::atan2(x, -1*y);            
			scalar xDir = Foam::cos(cellTheta);
			scalar yDir = Foam::sin(cellTheta);
			scalar r = Foam::sqrt(x*x + y*y)-R;

			Couette_[cellId] = (m*r+UMag)*vector(xDir,yDir,0);
		}
	}
	forAll(Couette_.mesh().boundary(), patchID)
	{
		forAll(Couette_.mesh().boundary()[patchID], faceId)
		{       
			if(unwrapped)
			{
				/* Unwrapped Journal Bearing */
				scalar y = Couette_.mesh().C().boundaryField()[patchID][faceId].y();
				Couette_.boundaryFieldRef()[patchID][faceId] = ((UMag/2)-m*y)*vector(1, 0, 0);
			}
			else
			{
				/* Journal Bearing Geometry */
				scalar y = Couette_.mesh().C().boundaryField()[patchID][faceId].y();
				scalar x = Couette_.mesh().C().boundaryField()[patchID][faceId].x();
				scalar faceTheta = Foam::atan2(x, -1*y);            
				scalar xDir = Foam::cos(faceTheta);
				scalar yDir = Foam::sin(faceTheta);
				scalar r = Foam::sqrt(x*x + y*y)-R;

				Couette_.boundaryFieldRef()[patchID][faceId] = (m*r+UMag)*vector(xDir, yDir, 0);
			}
		}
	}
}


// Calculates the film thickness based of the geometry of a bush bush
void TEHLThermo::calculateGeometric()
{
	scalar z0 = misalignmentCenter_.component(2).value();

	//Info << "Misalignment angle: " << misalignmentAngle_ << endl;
	//Info << "Misalignment center: " << misalignmentCenter_ << endl;

	scalar R = R_.value();

	// Loop through mesh and define geometric values for film thickness
	forAll(hGeo_, cellId)
	{
		scalar z 	= hGeo_.mesh().C()[cellId].component(2);
		scalar dz 	= z - z0;
		hGeo_[cellId] 	= c_.value() * (1-e*Foam::cos(theta_[cellId]+psi_)); // - dz*Foam::atan(misalignmentAngle_);   
	
	
		//scalar x	= hGeo_.mesh().C()[cellId].component(0);
		//hGeo_[cellId] 	= 1e-7 + R - Foam::sqrt(R*R - x*x);

		//hGeo_[cellId] = 2.54e-5+0.0176*Foam::pow(hGeo_.mesh().C()[cellId].x(), 2);
		//hGeo_[cellId] = 1.27e-5*(3-Foam::cos(82*hGeo_.mesh().C()[cellId].x()));
	}

	forAll(hGeo_.mesh().boundary(), patchID)
	{
		forAll(hGeo_.mesh().boundary()[patchID], faceId)
		{
			scalar z 	= hGeo_.mesh().C().boundaryField()[patchID][faceId].component(2);
			scalar dz 	= z - z0;

		    	hGeo_.boundaryFieldRef()[patchID][faceId] = c_.value() * (1-e*Foam::cos(theta_.boundaryFieldRef()[patchID][faceId]+psi_)); // - dz*Foam::atan(misalignmentAngle_); 
	
			//scalar x	= hGeo_.mesh().C().boundaryField()[patchID][faceId].component(0);
			//hGeo_.boundaryFieldRef()[patchID][faceId] = 1e-7 + 2*(R - Foam::sqrt(R*R - x*x));


			//hGeo_.boundaryFieldRef()[patchID][faceId] = 2.54e-5+0.0176*Foam::pow(hGeo_.mesh().C().boundaryField()[patchID][faceId].x(),2);
			//hGeo_.boundaryFieldRef()[patchID][faceId] = 1.27e-5*(3-Foam::cos(82*hGeo_.mesh().C().boundaryField()[patchID][faceId].x()));
		}
	}
};

// Calculates the a scalar field of the total surface deformation using the half-space assumption 
void TEHLThermo::calculateDeformation()
{
	if (Pstream::parRun())
	{
		List<vectorField> allCellCenters(Pstream::nProcs());
		List<scalarField> allTheta(Pstream::nProcs());
		List<scalarField> allCellVolumes(Pstream::nProcs());
		List<scalarField> allP(Pstream::nProcs());

		label myProc = Pstream::myProcNo();
		
		allCellCenters[myProc] = hDef_.mesh().C();	
		allCellVolumes[myProc] = hDef_.mesh().V();	
		allTheta[myProc] = theta_;
		allP[myProc] = p2D_.internalField();

		Pstream::gatherList(allCellCenters);
		Pstream::gatherList(allTheta);
		Pstream::gatherList(allCellVolumes);
		Pstream::gatherList(allP);
		Pstream::scatterList(allCellCenters);
		Pstream::scatterList(allTheta);
		Pstream::scatterList(allCellVolumes);
		Pstream::scatterList(allP);

		pressureIntegrationParallel(allP, allCellCenters, allCellVolumes, allTheta, hDef_, theta_, R_.value(), meshC_.value(), E_.value());
	}
	else
	{
		int nCells 	= hDef_.mesh().nCells();
		int cellStep 	= nCells/nThreads_;

		std::thread allThreads[int(nThreads_)];
		int cellStart, cellEnd;

		for (int i=0; i<nThreads_; i++)
		{
			cellStart = i*cellStep;

			if (i == (nThreads_-1))
			{
				cellEnd = nCells;
			}
			else
			{
				cellEnd = (i+1)*cellStep;
			}
			
			allThreads[i] = std::thread(pressureIntegration, 
							std::ref(p2D_),
							std::ref(hDef_),
							std::ref(theta_),
							R_.value(), 
							c_.value(), 
							E_.value(), 
							cellStart, 
							cellEnd);
					
		}

		for (int i=0; i<nThreads_; i++)
		{
			allThreads[i].join();
		}
		
		int cellId = 0;

		interpolateToBoundaries(hDef_);
	}
};

// Calcuated the thermal dilation using the temperature of the bush and bush surfaces
void TEHLThermo::calculateDilation()
{
	scalar shaftT(0);
	scalar bushT(0);
	scalar nFace(0);

	forAll(T_.mesh().boundary()[shaftPatchId_], faceId)
	{
		nFace += 1;
		shaftT += T2D_.boundaryFieldRef()[shaftPatchId_][faceId];
	}

	forAll(T_.mesh().boundary()[bushPatchId_], faceId)
	{
		bushT += T2D_.boundaryFieldRef()[bushPatchId_][faceId];
	}

	if (Pstream::parRun())
	{
		reduce(shaftT, sumOp<scalar>());
		reduce(bushT, sumOp<scalar>());		
		reduce(nFace, sumOp<scalar>());
	}

	shaftT 	= shaftT/nFace;
	bushT 	= bushT/nFace;

	hDil_.value() 	= alphaShaft_.value()*(R_.value()+c_.value())*(shaftT - TRef_.value()) - alphaBush_.value()*R_.value()*(bushT - TRef_.value());

	Info << "Journal dilation = " << hDil_.value() *1e6 << " um" << endl;
}

// Calculates a viscosity field based on the temperature of the fluid
void TEHLThermo::calculateMu()
{
	dimensionedScalar T1(dimless/dimTemperature, muTCoeff_[0].value());
	dimensionedScalar T2(dimTemperature, muTCoeff_[1].value());
	
	scalar P1 = muPCoeff_[0].value();
	scalar P2 = muPCoeff_[1].value();
	scalar P3 = muPCoeff_[2].value();
	
	scalar muOil = muConst_.value();
	scalar muAir = muAirConst_.value();

	volScalarField mu_l = muConst_*Foam::exp(T1*(T2D_ - T2));
	volScalarField mu_eqv = (1-relRho2D_)*muAirConst_ + relRho2D_*mu_l;

	mu2D_ = (1-g_)*mu_eqv + g_*mu_l;
	
	/*
	forAll(mu2D_, cellId)
	{
		scalar ff = relRho2D_[cellId];

		// Exponential viscosity function of temperature with cavitation effects 
		if (ff > 1)
		{
			ff = 1;
		}

		scalar muT = muOil*Foam::exp(T1*(T2D_[cellId] - T2));

		//Exponential viscosity function including pressure and cavitation effects 

		//mu2D_[cellId] = muAir + (muT*Foam::exp(P1*Foam::exp(P2/((T2D_[cellId]-273)+P3))*(p2D_[cellId]+1e5)*1e-9) - muAir)*ff;
	}	

	forAll(mu2D_.mesh().boundary(), patchID)
	{
		forAll(mu2D_.mesh().boundary()[patchID], faceId)
		{
			scalar ff = relRho2D_.boundaryFieldRef()[patchID][faceId];

			if (ff > 1)
			{
				ff = 1;
			}

			//mu_.boundaryFieldRef()[patchID][faceId] = muAir + (muOil*Foam::exp(T1*(T_.boundaryFieldRef()[patchID][faceId] - T2)) - muAir)*filmFraction_.boundaryFieldRef()[patchID][faceId];
			//mu2D_.boundaryFieldRef()[patchID][faceId] = muOil*Foam::exp(T1*(T2D_.boundaryFieldRef()[patchID][faceId] - T2));

			scalar muT = muOil*Foam::exp(T1*(T2D_.boundaryFieldRef()[patchID][faceId] - T2));

			mu2D_.boundaryFieldRef()[patchID][faceId] = muAir + (muT*Foam::exp(P1*Foam::exp(P2/((T2D_.boundaryFieldRef()[patchID][faceId]-273)+P3))*(p2D_.boundaryFieldRef()[patchID][faceId]+1e5)*1e-9) - muAir)*ff;

		}
	}
	*/
};

// Calculates a field on Bulk Modulus values
void TEHLThermo::calculateB()
{
	scalar dPhi = 0.05;
	scalar BGrad = (BConst_.value()-BCav_.value())/dPhi;

	/**** Linear Variation for Conventional Elrod-Adams Model *****/
	forAll(relRho2D_, cellId)
	{
		if (relRho2D_[cellId] >= (1+dPhi))
		{
			B_[cellId] = BConst_.value();
		}
		else if (relRho2D_[cellId] < 1)
		{
			B_[cellId] = BCav_.value();
		}
		else
		{
			B_[cellId] = BCav_.value() + BGrad*(relRho2D_[cellId] - 1); 
		}
	}
};

// Calculates pressure for different cavitation models
void TEHLThermo::calculateP()
{
	// Calculate pressure from relRho for Elrod Adams Model
	forAll(relRho2D_, cellId)
	{
		if (relRho2D_[cellId] < 1)
		{
			p2D_[cellId] = PCav_.value();
		}
		else
		{
			p2D_[cellId] = PCav_.value() + B_[cellId]*Foam::log(relRho2D_[cellId]);
		}
	}

	forAll(relRho2D_.mesh().boundary(), patchId)
	{
		const scalarField& relRhoPatch = relRho2D_.boundaryField()[patchId];
		scalarField& pPatch = p2D_.boundaryFieldRef()[patchId];

		forAll(relRho2D_.mesh().boundary()[patchId], faceId)
		{
			if (relRhoPatch[faceId] < 1)
			{
				pPatch[faceId] = PCav_.value();
			}
			else
			{
				pPatch[faceId] = PCav_.value() + BConst_.value()*Foam::log(relRhoPatch[faceId]);
			}
		}
	}

};

// Applies the Elrod Adam switch function
void TEHLThermo::elrodAdamSwitch()
{	
 	/******* Switch Function incorporating modification from Fesanghary et al *******/
	forAll(g_, cellId)
	{
		if (relRho2D_[cellId] >= 1)
		{
			if (gFactor_ > 0)
			{
			    	g_[cellId] = (g_[cellId]/gFactor_)*relRho2D_[cellId];
				//g_[cellId] = (g_[cellId]/gFactor_);
			}
			else
			{
			    g_[cellId] = 1;
			}
		}
		else
		{
			g_[cellId] = g_[cellId] * gFactor_;
		}

		if (g_[cellId] > 1)
		{
			g_[cellId] = 1;
		}
		if (g_[cellId] < gLimit_)
		{
			g_[cellId] = gLimit_;
		}
		/*
		if (g_[cellId] < 1e-6)
		{
			g_[cellId] = 0;
		}
		*/
	}
};

void TEHLThermo::setSwitchFunction(scalar value)
{
	if(value < 0)
	{
		value = gLimit_;
	}

	forAll(g_, cellId)
	{
		if (g_[cellId] < 1)
		{
			g_[cellId] = value;
		}
		else
		{
			g_[cellId] = 1;
		}
	}
};

// Applies different cavitation methods
void TEHLThermo::calculateCavFields()
{
	/**** Update fields for Elrod-Adam *****/
	elrodAdamSwitch();	
	calculateP();
	//calculateB();
};

// Cacluates new film thickness field
void TEHLThermo::updateFilmThickness()
{
	h2D_.storePrevIter();

	if(deformation)
	{
		calculateDeformation();
	}
	if(dilation)
	{
		calculateDilation();
	}

	h2D_ = hGeo_ + hDef_ + hDil_;
	h2D_.relax();
};

// Resets empirically calculated fields for restarting the numerical procedure
void TEHLThermo::reset()
{
	calculateGeometric();
	h2D_ = hGeo_ + hDef_ + hDil_;
	
	forAll(g_, cellId)
	{
		g_[cellId] = 1;
	}
};

// Calculate the Pouseille flow field and update velocity field USum
void TEHLThermo::updateVelocity()
{
	volVectorField gradPhi = fvc::grad(relRho_);
	volVectorField gradP = fvc::grad(p_);

	scalar c = meshC_.value();
	scalar R = R_.value();
	scalar B = BConst_.value();

	
	if(includePoiseuille)
	{
		/* Calculate the Poiseuille velocity profile from the pressure solution */
		forAll(Poiseuille_, cellId)
		{  
			if(unwrapped)
			{
				scalar y 	= Poiseuille_.mesh().C()[cellId].y()+(c/2);
				scalar g 	= g3D_[cellId];
				scalar mu 	= mu_[cellId];
				scalar h 	= h_[cellId];
				scalar hc 	= h/c;

				Poiseuille_[cellId] = (Foam::pow(hc,2)*y*y - h*hc*y)*gradP[cellId]/(2*mu);
			}
			else
			{
				scalar cellx 	= Poiseuille_.mesh().C()[cellId].x();
				scalar celly 	= Poiseuille_.mesh().C()[cellId].y();
				scalar y 	= Foam::sqrt(cellx*cellx+celly*celly) - R;
				scalar g 	= g3D_[cellId];		
				scalar mu 	= mu_[cellId];
				scalar h 	= h_[cellId];
				scalar hc 	= h/c;

				Poiseuille_[cellId] = (Foam::pow(hc,2)*y*y - h*hc*y)*gradP[cellId]/(2*mu);
			}
		}

		/* Calculate for the boundary faces */
		forAll(Poiseuille_.mesh().boundary(), patchID)
		{
			forAll(Poiseuille_.mesh().boundary()[patchID], cellId)
			{       
				if(unwrapped)
				{
					scalar y 	= Poiseuille_.mesh().C().boundaryField()[patchID][cellId].y()+(c/2);
					scalar g 	= g3D_.boundaryFieldRef()[patchID][cellId];
					scalar mu 	= mu_.boundaryFieldRef()[patchID][cellId];
					scalar h 	= h_.boundaryFieldRef()[patchID][cellId];
					scalar hc 	= h/c;

					Poiseuille_.boundaryFieldRef()[patchID][cellId] =  (Foam::pow(hc,2)*y*y - h*hc*y)*gradP.boundaryFieldRef()[patchID][cellId]/(2*mu);	
				}
				else
				{

					scalar cellx 	= Poiseuille_.mesh().C().boundaryField()[patchID][cellId].x();
					scalar celly 	= Poiseuille_.mesh().C().boundaryField()[patchID][cellId].y();
					scalar y 	= Foam::sqrt(cellx*cellx+celly*celly) - R;
					scalar g 	= g3D_.boundaryFieldRef()[patchID][cellId];
					scalar mu 	= mu_.boundaryFieldRef()[patchID][cellId];
					scalar h 	= h_.boundaryFieldRef()[patchID][cellId];
					scalar hc 	= h/c;

					Poiseuille_.boundaryFieldRef()[patchID][cellId] =  (Foam::pow(hc,2)*y*y - h*hc*y)*gradP.boundaryFieldRef()[patchID][cellId]/(2*mu);	
				}
			}
		}

		/* Full velocity profile calcualte as sum of Couette and Poiseuille fields */
		U_ = (Couette_ + Poiseuille_);
	}
	else
	{
		/* Excludes Poiseuille field, for 2D cases */
		U_ = Couette_;
	}
};

/* Calculate load in the vertical direction */
dimensionedScalar TEHLThermo::calculateLoadVertical()
{
	scalar Area;
	vector force(0,0,0);
	dimensionedScalar load(dimensionSet(1,1,-2,0,0,0,0), 0);

	forAll(p_, cellID)
	{
		Area = p_.mesh().V()[cellID]/c_.value(); 
		force += (p_[cellID]-PCav_.value())*Area*vector(-1*Foam::sin(theta_[cellID]), Foam::cos(theta_[cellID]), 0);
	}

	// Force assumed in vertical
	load.value() = (force & vector(0,1,0));

	return load;
};

// Calculate load and angular offset of the bush
dimensionedScalar TEHLThermo::calculateLoad()
{
	vector force(0,0,0);
	dimensionedScalar load(dimensionSet(1,1,-2,0,0,0,0), 0);


	if(unwrapped)
	{
		const labelUList& cellList = mesh2D_.boundaryMesh()[shaftPatchId_].faceCells();
		const scalarField& faceAreas = mesh2D_.magSf().boundaryField()[shaftPatchId_];

		forAll(cellList, faceId)
		{
			scalar& ang = theta_[cellList[faceId]];
			scalar& pRef = PCav_.value();
			vector cellMagVec = vector(Foam::sin(ang), Foam::cos(ang), 0) * faceAreas[faceId];

			force += (p2D_[cellList[faceId]]-pRef) * cellMagVec;
		}
	}
	else
	{
		const labelUList& cellList = mesh2D_.boundaryMesh()[shaftPatchId_].faceCells();
		vectorField faceVectors = mesh2D_.Sf().boundaryField()[shaftPatchId_];

		forAll(faceVectors, faceId)
		{
			force += (p2D_[cellList[faceId]]-PCav_.value())*faceVectors[faceId];
		}
	}

	// Load assumed vertical and locked angular offset
	//load.value() = (force & vector(0,1,0));
		
	// Varies angular offset
	loadAngle_ = Foam::acos((force & vector(0, 1, 0))/mag(force));
	
	if (force.x() > 0)
	{
		psi_ -= loadAngle_;
	}
	else
	{
		psi_ += loadAngle_;
	}

	if (force.component(0) < 0)
	{
		loadAngle_ = -1*loadAngle_;
	}

	load.value() = mag(force);
	
	return load;
};

// Calculate load and angular offset of the bush
vector TEHLThermo::calculateForce()
{
	vector force(0,0,0);

	const fvMesh& mesh = p2D_.mesh();

	const labelUList& cellList = mesh2D_.boundaryMesh()[shaftPatchId_].faceCells();
	vectorField faceVectors = mesh2D_.Sf().boundaryField()[shaftPatchId_];
	scalarField pShaft = p2D_.boundaryFieldRef()[shaftPatchId_];

	forAll(pShaft, faceId)
	{
		//force += (p2D_[cellList[faceId]]-PCav_.value())*faceVectors[faceId];
		force += pShaft[faceId] * faceVectors[faceId]; 
	}

	return force;
};

void TEHLThermo::updatePsi(scalar dPsi)
{
	psi_ += dPsi;
};

void TEHLThermo::calculateAddressing()
{
	forAll(mesh3D_.C(), cellId)
	{
		const vector& cellCenter = mesh3D_.C()[cellId];
		label nearestCellId = mesh2D_.findNearestCell(cellCenter);
		cellAddressing_[nearestCellId].append(cellId);
	}

	forAll(mesh2D_.Cf().boundaryField(), patchId)
	{
		labelListList boundaryAddressing(mesh2D_.Cf().boundaryField()[patchId].size());

		const vectorField& boundaryCenters = mesh2D_.Cf().boundaryField()[patchId];
		
		forAll(mesh3D_.boundaryMesh()[patchId], faceId)
		{
			const vector& faceCenter = mesh3D_.Cf().boundaryField()[patchId][faceId];

			label nearestFaceId = nearestBoundaryFace(faceCenter, boundaryCenters);
			
			boundaryAddressing[nearestFaceId].append(faceId);		
		}

		faceAddressing_[patchId].append(boundaryAddressing);
	}
};

void TEHLThermo::interpolateT()
{
	Info << "Interpolating T" << endl;
	forAll(mesh2D_.C(), cellId)
	{
		List<label> cellList = cellAddressing_[cellId];
		scalar TSum(0);
		scalar nCells(0);
		scalar volSum(0);

		forAll(cellList, i)
		{
			scalar vol = mesh3D_.V()[cellList[i]];
			TSum += T_[cellList[i]]*vol;
			volSum += vol;
			nCells ++;
		}

		T2D_[cellId] = TSum/volSum;
	}		

	interpolateToBoundaries(T2D_);
};

void TEHLThermo::interpolateTo3D()
{
	forAll(T2D_, cellId)
	{
		List<label> cellList = cellAddressing_[cellId];

		forAll(cellList, i)
		{
			relRho_[cellList[i]] = relRho2D_[cellId];
			p_[cellList[i]] = p2D_[cellId];
			g3D_[cellList[i]] = g_[cellId];
			h_[cellList[i]] = h2D_[cellId];
			mu_[cellList[i]] = mu2D_[cellId];
		}
	}		

	forAll(T2D_.boundaryField(), patchId)
	{
		scalarField& relRhoBoundary = relRho2D_.boundaryFieldRef()[patchId];
		scalarField& pBoundary = p2D_.boundaryFieldRef()[patchId];

		forAll(relRhoBoundary, faceId)
		{
			List<label> faceList = faceAddressing_[patchId][faceId];

			forAll(faceList, i)
			{
				relRho_.boundaryFieldRef()[patchId][faceList[i]] = relRho2D_.boundaryField()[patchId][faceId];
				p_.boundaryFieldRef()[patchId][faceList[i]] = p2D_.boundaryField()[patchId][faceId];
			}
		}
	}

	interpolateToBoundaries(g_);
	interpolateToBoundaries(h_);
	interpolateToBoundaries(mu_);
};

void TEHLThermo::interpolatePToRelRho()
{
	scalar pCav = PCav_.value();
	scalar B = BConst_.value();

	forAll(interpBoundaries_, boundaryNameId)
	{
		word& boundaryName = interpBoundaries_[boundaryNameId];
		label patchId2D = mesh2D_.boundaryMesh().findPatchID(boundaryName);
		label patchId3D = mesh3D_.boundaryMesh().findPatchID(boundaryName);
		
		forAll(mesh2D_.boundaryMesh()[patchId2D], faceId)
		{
			List<label> faceList = faceAddressing_[patchId2D][faceId];

			scalar nCells(0);
			scalar pSum(0);

			forAll(faceList, i)
			{
				nCells ++;
				pSum += p_.boundaryField()[patchId3D][faceList[i]];
			}

			pSum = pSum/nCells;

			scalar faceRelRho = Foam::exp((pSum-pCav)/B);

			if (faceRelRho < 1)
			{
				relRho_.boundaryFieldRef()[patchId2D][faceId] = 1;
			}
			else
			{
				relRho_.boundaryFieldRef()[patchId2D][faceId] = faceRelRho;
			}
		}
	}
};

void TEHLThermo::viscosityIntergals()
{
	scalar R = R_.value();

	forAll(F2_, cellId)
	{
		List<label> cellList = cellAddressing_[cellId];
		scalar F0(0);
		scalar F1(0);
		scalar F2(0);

		scalar dy = meshC_.value()/cellList.size();

		// Calculate F0 and F1
		forAll(cellList, i)
		{
			scalar mu = mu_[cellList[i]];
			F0 += dy/mu;
			scalar y = mu_.mesh().C()[cellList[i]].y();
			scalar x = mu_.mesh().C()[cellList[i]].x();
			scalar r = Foam::sqrt(x*x + y*y)-R;

			F1 += dy*r/mu;
		}

		scalar F1F0 = F1/F0;

		// Calculate F2
		forAll(cellList, i)
		{
			scalar mu = mu_[cellList[i]];
			scalar y = mu_.mesh().C()[cellList[i]].y();
			scalar x = mu_.mesh().C()[cellList[i]].x();
			scalar r = Foam::sqrt(x*x + y*y)-R;
			
			F2 += dy*((r/mu)*(r - F1F0));
		}

		F2_[cellId] = F2;
		F3_[cellId] = h2D_[cellId]*(1-F1F0);
	}	

	interpolateToBoundaries(F2_);
	interpolateToBoundaries(F3_);
};

void TEHLThermo::updateRho()
{
	volScalarField rho_eqv = (1 - relRho_)*rhoAirConst_ + relRho_*rhoConst_;

	rho_ = (1-g3D_)*rho_eqv + g3D_*rhoConst_;
	/*
	forAll(rho_, cellId)
	{
		scalar& relRho = relRho_[cellId];
		scalar& rhof = rhoConst_.value();
		scalar& rhog = rhoAirConst_.value();

		if (relRho > 1)
		{
			rho_[cellId] = relRho*rhof;
		}
		else
		{
			rho_[cellId] = relRho*rhof + (1-relRho)*rhog;
		}
	}
	*/
};	

// Calculate virtual thermal conductivity to account for the difference between film thickness and the mesh geometry
void TEHLThermo::updatek()
{
	//k_ = kAirConst_ + (kConst_ - kAirConst_)*filmFraction_;
	k_ = 1/((filmFraction_/kConst_) + ((1-filmFraction_)/kAirConst_));
	volScalarField k_eqv = (1-relRho_)*kAirConst_ + relRho_*kConst_;
	k_ = (1 - g3D_)*k_eqv + g3D_*kConst_;
};	

void TEHLThermo::updatecp()
{
	volScalarField cp_eqv = ((1-relRho_)*rhoAirConst_*cpAirConst_ + relRho_*rhoConst_*cpConst_)/((1-relRho_)*rhoAirConst_ + relRho_*rhoConst_);
	cp_ = (1-g3D_)*cp_eqv + g3D_*cpConst_;	
	//cp_ = cpAirConst_ + (cpConst_ - cpAirConst_)*filmFraction_;
};	

void TEHLThermo::updateFilmFraction()
{
	forAll(filmFraction_, cellId)
	{
		if (relRho_[cellId] < 1)
		{
			filmFraction_[cellId] = relRho_[cellId];
		}
		else
		{
			filmFraction_[cellId] = 1;
		}
	}

	forAll(filmFraction_.mesh().boundary(), patchID)
	{
		forAll(filmFraction_.mesh().boundary()[patchID], faceId)
		{
			if (relRho_.boundaryFieldRef()[patchID][faceId] < 1)
			{
				filmFraction_.boundaryFieldRef()[patchID][faceId] = relRho_.boundaryFieldRef()[patchID][faceId];
			}
			else
			{
				filmFraction_.boundaryFieldRef()[patchID][faceId] = 1;
			}
		}
	}
};	

/************************************* Access Functions****************************************/

volScalarField& TEHLThermo::p()
{
	return p_;
};

volScalarField& TEHLThermo::p2D()
{
	return p2D_;
};

volScalarField& TEHLThermo::T()
{
	return T_;
};

volScalarField& TEHLThermo::T2D()
{
	return T2D_;
};

volScalarField& TEHLThermo::relRho()
{
	return relRho_;
};

volScalarField& TEHLThermo::relRho2D()
{
	return relRho2D_;
};

volScalarField& TEHLThermo::filmFraction()
{
	return filmFraction_;
};

volScalarField& TEHLThermo::g()
{
	return g_;
};

volScalarField& TEHLThermo::g3D()
{
	return g3D_;
};

volScalarField& TEHLThermo::B()
{
	return B_;
};

volScalarField& TEHLThermo::h()
{
	return h_;
};

volScalarField& TEHLThermo::h2D()
{
	return h2D_;
};

volScalarField& TEHLThermo::hDef()
{
	return hDef_;
};

volScalarField& TEHLThermo::hGeo()
{
	return hGeo_;
};

volScalarField& TEHLThermo::theta()
{
	return theta_;
};

volScalarField& TEHLThermo::rhoFilm()
{
	return rho_;
};

volScalarField& TEHLThermo::muFilm()
{
	return mu_;
};

volScalarField& TEHLThermo::muFilm2D()
{
	return mu2D_;
};

volVectorField& TEHLThermo::radial()
{
	return radial_;
};

volVectorField& TEHLThermo::Couette()
{
	return Couette_;
};

volVectorField& TEHLThermo::Poiseuille()
{
	return Poiseuille_;
};

volVectorField& TEHLThermo::U()
{
	return U_;
};

volVectorField& TEHLThermo::U2D()
{
	return U2D_;
};

volScalarField TEHLThermo::hByc()
{
	volScalarField hByc
	(
		IOobject
		(
			"hByc",
			h_.mesh().time().timeName(),
			h_.mesh(),
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		h_/meshC_
	);

	return hByc;
};

dimensionedScalar& TEHLThermo::PCav()
{
	return PCav_;
};

dimensionedScalar& TEHLThermo::c()
{
	return c_;
};

dimensionedScalar& TEHLThermo::R()
{
	return R_;
};

dimensionedScalar& TEHLThermo::L()
{
	return L_;
};

dimensionedScalar& TEHLThermo::meshC()
{
	return meshC_;
};

dimensionedScalar& TEHLThermo::dRhoT()
{
	return dRhoT_;
};

volScalarField& TEHLThermo::cp()
{
	return cp_;
};

volScalarField& TEHLThermo::F2()
{
	return F2_;
};

volScalarField& TEHLThermo::F3()
{
	return F3_;
};

surfaceScalarField TEHLThermo::massFlux()
{
	surfaceScalarField massFlux_ = fvc::flux(relRho_*U_);

	return massFlux_;
};

surfaceScalarField TEHLThermo::kSurf()
{
	this->updatek();
	surfaceScalarField kSurf_ = mag(fvc::flux(k_*radial_/this->hByc())/mesh3D_.magSf());

	return kSurf_;
};

dimensionedScalar& TEHLThermo::rhoConst()
{
	return rhoConst_;
};

dimensionedScalar& TEHLThermo::muConst()
{
	return muConst_;
};

dimensionedScalar& TEHLThermo::omega()
{
	return omega_;
};

dimensionedScalar TEHLThermo::UShaft()
{
	return omega_*R_;
};

dimensionedScalar& TEHLThermo::bushT()
{
	return bushT_;
};

dimensionedScalar& TEHLThermo::hDil()
{
	return hDil_;
};

scalar& TEHLThermo::psi()
{
	return psi_;
};

scalar& TEHLThermo::loadAngle()
{
	return loadAngle_;
};
   	    	
bool& TEHLThermo::includeDeformation()
{
	return deformation;
};

bool& TEHLThermo::includeTemperature()
{
	return temperature;
};

bool& TEHLThermo::includeMisalignment()
{
	return misalignment;
};

bool& TEHLThermo::unwrappedCase()
{
	return unwrapped;
};

label& TEHLThermo::shaftPatchId()
{
	return shaftPatchId_;
};

label& TEHLThermo::bushPatchId()
{
	return bushPatchId_;
};

label& TEHLThermo::sidesPatchId()
{
	return sidesPatchId_;
};

label& TEHLThermo::inletPatchId()
{
	return inletPatchId_;
};

label& TEHLThermo::outletPatchId()
{
	return outletPatchId_;
};

/**************** Derived Functions ****************/
void TEHLThermo::correct(){};

word TEHLThermo::thermoName() const 
{
	return phaseName_;
};	

bool TEHLThermo::incompressible() const 
{
	return false;
};

bool TEHLThermo::isochoric() const 
{
	return false;
};

tmp<volScalarField> TEHLThermo::rho() const
{
	return p_;
};

tmp<scalarField> TEHLThermo::rho(const label patchi) const
{
	return p_.boundaryField()[patchi];
};

volScalarField& TEHLThermo::he()
{
	return p_;
};

tmp<volScalarField> TEHLThermo::he(const volScalarField& p, const volScalarField& T) const
{	
	return p_;
};

tmp<scalarField> TEHLThermo::he(const scalarField& T,const labelList& cells) const
{	
	return T;
};

tmp<scalarField> TEHLThermo::he(const scalarField& T,const label patchi) const
{
	return p_.boundaryField()[patchi];
};

const volScalarField& TEHLThermo::he() const 
{
	return p_;
};

tmp<volScalarField> TEHLThermo::hs() const
{
	return p_;
};

tmp<volScalarField> TEHLThermo::hs(const volScalarField& p,const volScalarField& T) const
{
	return p_;
};

tmp<scalarField> TEHLThermo::hs(const scalarField& T,const labelList& cells) const
{	
	return T;
};

tmp<scalarField> TEHLThermo::hs(const scalarField& T,const label patchi) const
{
	return p_.boundaryField()[patchi];
};

tmp<volScalarField> TEHLThermo::ha() const 
{
	return p_;
};

tmp<volScalarField> TEHLThermo::ha(const volScalarField& p,const volScalarField& T) const
{
	return p_;
};

tmp<scalarField> TEHLThermo::ha(const scalarField& T,const labelList& cells) const
{	
	return T;
};

tmp<scalarField> TEHLThermo::ha(const scalarField& T,const label patchi) const
{
	return p_.boundaryField()[patchi];
};

tmp<scalarField> TEHLThermo::THE(const scalarField& h,const scalarField& T0,const labelList& cells) const
{
	return h;
};

tmp<scalarField> TEHLThermo::THE(const scalarField& h, const scalarField& T0,const label patchi) const
{
	return p_.boundaryField()[patchi];
};

tmp<volScalarField> TEHLThermo::hc() const 
{
	return p_;
};

tmp<volScalarField> TEHLThermo::Cp() const 
{
	return p_;
};

tmp<scalarField> TEHLThermo::Cp(const scalarField& T,const label patchi) const
{
	return p_.boundaryField()[patchi];
};

tmp<volScalarField> TEHLThermo::Cv() const 
{
	return p_;
};

tmp<scalarField> TEHLThermo::Cv(const scalarField& T,const label patchi) const
{
	return p_.boundaryField()[patchi];	
};

tmp<volScalarField> TEHLThermo::Cpv() const 
{
	return p_;
};

tmp<scalarField> TEHLThermo::Cpv(const scalarField& T,const label patchi) const
{
	return p_.boundaryField()[patchi];
};

tmp<volScalarField> TEHLThermo::CpByCpv() const 
{
	return p_;
};

tmp<scalarField> TEHLThermo::CpByCpv(const scalarField& T,const label patchi) const
{
	return p_.boundaryField()[patchi];
};

tmp<volScalarField> TEHLThermo::kappa() const
{
	return k_;
};

tmp<scalarField> TEHLThermo::kappa(const label patchi) const
{
	return k_.boundaryField()[patchi];
};

tmp<volScalarField> TEHLThermo::alphahe() const 
{
	return p_;
};

tmp<scalarField> TEHLThermo::alphahe(const label patchi) const 
{
	return p_.boundaryField()[patchi];
};

void TEHLThermo::correctRho(const volScalarField& deltaRho){};

const volScalarField& TEHLThermo::psi() const 
{
	return p_;
};

tmp<volScalarField> TEHLThermo::gamma() const 
{
	return p_;
};

tmp<scalarField> TEHLThermo::gamma(const scalarField& T,const label patchi) const
{
	return p_.boundaryField()[patchi];
};

tmp<volScalarField> TEHLThermo::W() const 
{
	return p_;
};

tmp<scalarField> TEHLThermo::W(const label patchi) const 
{
	return p_.boundaryField()[patchi];
};

tmp<volScalarField> TEHLThermo::mu() const 
{
	return p_;
};

tmp<scalarField> TEHLThermo::mu(const label patchi) const 
{
	return p_.boundaryField()[patchi];
};

tmp<volScalarField> TEHLThermo::kappaEff(const volScalarField&) const
{
	return p_;	
};

tmp<scalarField> TEHLThermo::kappaEff(const scalarField& alphat,const label patchi) const
{
	return p_.boundaryField()[patchi];	
};

tmp<volScalarField> TEHLThermo::alphaEff(const volScalarField& alphat) const
{
	return p_;	
};

tmp<scalarField> TEHLThermo::alphaEff(const scalarField& alphat,const label patchi) const
{
	return p_.boundaryField()[patchi];
};

