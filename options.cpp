#include "options.h"




options::options(int argc, char** argv) :input(""), output(""),
winSize(75), gapColFrac(1.f),
avgToOthers(false), avgToOthersFile(""),fullDistMat(false), distMatOut(""),
minGoodPosFrac(0.f),
NT(true),AA(false),quiet(true),
maskBorderGap(false), maskLowID(false) {
	if (argc == 0) { return; }//could be the pseudo mod opt object

	for (int i = 1; i < argc; i++)
	{
		if (!strcmp(argv[i], "-i"))
			input = argv[++i];
		else if (!strcmp(argv[i], "-o"))
			output = argv[++i];
		else if (!strcmp(argv[i], "-avgID"))   //get the average ID to all other seqs to target
			avgToOthersFile = argv[++i];
		else if (!strcmp(argv[i], "-noFullDist"))   //calculate full distance matrix, file to save this
			distMatOut = argv[++i];
		else if (!strcmp(argv[i], "-nt"))   
			NT = true;
		else if (!strcmp(argv[i], "-aa"))   
			AA = false;
		else if (!strcmp(argv[i], "-gapsNoPen"))   
			GapsPenalty = false;
		else if (!strcmp(argv[i], "-maskBorderGap"))  
			maskBorderGap = true;
		else if (!strcmp(argv[i], "-maskLowID"))   
			maskLowID = true;
		else if (!strcmp(argv[i], "-rmGapColsGreater"))
			gapColFrac = atof(argv[++i]);
		else if (!strcmp(argv[i], "-minGoodPosFrac"))
			minGoodPosFrac = atof(argv[++i]);
		

	}
	if (avgToOthersFile != "" || maskLowID) {
		avgToOthers = true;
	}
	if (distMatOut != "" || avgToOthers) {
		fullDistMat = true;
	}

	if (AA) { NT = false; }
	if (input == output || output == "") { 
		cerr << "Overwriting input file\n"; 
		output = input;
	}

}