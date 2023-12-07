class inputTEHL
{
private:

	const IOdictionary settingsDict;
	const IOdictionary filmControlDict;
	const IOdictionary solidConstantsDict;
	const IOdictionary fluidConstantsDict;
	const IOdictionary outputControlDict;

public:

	// Constructor
	inputTEHL
	(
		fvMesh& filmMesh,
		Time& runTime
	);
    
    	// Destructor
    	~inputTEHL();
    	
    	IOdictionary settings();
    	
	IOdictionary filmControl();
	
	IOdictionary solidConstants();
	
	IOdictionary fluidConstants();
	
	IOdictionary outputControl();

};
