class controlOutputTEHL
{
private:

	bool relRho;
	bool filmFraction;
	bool g;
	bool B;
	bool h;
	bool hRigid;
	bool hDeformation;
	bool hDilation;
	bool rho;
	bool mu;
	bool couette;
	bool poiseuille;
	bool U;
	bool k;
	bool cp;

	TEHLThermo* filmThermo;
	Time* runTime;

public:

	// Constructor
	controlOutputTEHL
	(
		IOdictionary controlDict,
		TEHLThermo& thermo,
		Time& time
	);
    
    	// Destructor
    	~controlOutputTEHL();

    	// Write all variables requested
    	void write();

	void writeNow();
};
