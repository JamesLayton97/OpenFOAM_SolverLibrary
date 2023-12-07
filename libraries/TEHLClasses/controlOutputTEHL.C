#include "TEHLClasses.H"

controlOutputTEHL::controlOutputTEHL
(
	IOdictionary controlDict,
	TEHLThermo& thermo,
	Time& time
)
{

	relRho 			= readBool(controlDict.lookup("relRho"));
	filmFraction 		= readBool(controlDict.lookup("filmFraction"));
	g 			= readBool(controlDict.lookup("g"));
	B 			= readBool(controlDict.lookup("bulkModulus"));
	h			= readBool(controlDict.lookup("filmThickness"));
	hRigid 			= readBool(controlDict.lookup("rigidFilmThickness"));
	hDeformation 		= readBool(controlDict.lookup("deformation"));
	rho 			= readBool(controlDict.lookup("rho"));
	mu 			= readBool(controlDict.lookup("mu"));
	couette 		= readBool(controlDict.lookup("couette"));
	poiseuille 		= readBool(controlDict.lookup("poiseuille"));
	U			= readBool(controlDict.lookup("U"));
	k 			= readBool(controlDict.lookup("k"));
	cp 			= readBool(controlDict.lookup("cp"));

	filmThermo = &thermo;
	runTime = &time;
};

controlOutputTEHL::~controlOutputTEHL()
{};

void controlOutputTEHL::write()
{
	if (runTime->write())
	{
		// Write fields if defined true
		if (relRho)
		{
			filmThermo->relRho2D().write();
		}
		if (filmFraction)
		{
			filmThermo->filmFraction().write();
		}
		if (g)
		{
			filmThermo->g().write();
		}
		if (B)
		{
			filmThermo->B().write();
		}
		if (h)
		{
			filmThermo->h2D().write();
		}
		if (hRigid)
		{
			filmThermo->hGeo().write();
		}
		if (hDeformation)
		{
			filmThermo->hDef().write();
		}
		if (rho)
		{
			filmThermo->rhoFilm().write();
		}
		if (mu)
		{
			filmThermo->muFilm2D().write();
			filmThermo->muFilm().write();
		}
		if (couette)
		{
			filmThermo->Couette().write();
		}
		if (poiseuille)
		{
			filmThermo->Poiseuille().write();
		}
		if (U)
		{
			filmThermo->U().write();  // Error causing divergence in T solution
		}
		if (k)
		{
			volScalarField k = filmThermo->kappa();
			k.write();
		}
		if (cp)
		{
			filmThermo->cp().write();
		}
	}
};		

void controlOutputTEHL::writeNow()
{
	// Write fields if defined true
	if (relRho)
	{
		filmThermo->relRho2D().write();
	}
	if (filmFraction)
	{
		filmThermo->filmFraction().write();
	}
	if (g)
	{
		filmThermo->g().write();
	}
	if (B)
	{
		filmThermo->B().write();
	}
	if (h)
	{
		filmThermo->h2D().write();
	}
	if (hRigid)
	{
		filmThermo->hGeo().write();
	}
	if (hDeformation)
	{
		filmThermo->hDef().write();
	}
	if (rho)
	{
		filmThermo->rhoFilm().write();
	}
	if (mu)
	{
		filmThermo->muFilm2D().write();
		filmThermo->muFilm().write();
	}
	if (couette)
	{
		filmThermo->Couette().write();
	}
	if (poiseuille)
	{
		filmThermo->Poiseuille().write();
	}
	if (U)
	{
		filmThermo->U().write();  // Error causing divergence in T solution
	}
	if (k)
	{
		volScalarField k = filmThermo->kappa();
		k.write();
	}
	if (cp)
	{
		filmThermo->cp().write();
	}
};		
