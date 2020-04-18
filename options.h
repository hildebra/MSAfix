#pragma once
#include <stdio.h>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <assert.h>
#include <algorithm>
#include "string.h"

using namespace std;
typedef  unsigned int uint;


struct options
{
public:
	options(int argc, char** argv);
	//vars
	std::string input = "";
	std::string output = "";
	uint winSize = 150;
	float gapColFrac;
	bool avgToOthers = false; //get the average ID to all other seqs to target
	string avgToOthersFile;
	bool fullDistMat = true;
	string distMatOut;
	float minGoodPosFrac;
	bool GapsPenalty = true;
	bool NT, AA;
	bool quiet;
	bool maskBorderGap;
	bool maskLowID;
};