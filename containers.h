#pragma once
#include "options.h"

//DNA constants (for comparisons)
static const char DNA_SPACE[15] = { 'A','C','G','T','N','R','Y','M','K','W','S','B','D','H','V' };
static const char AA_SPACE[25] = { 'A','B','C','D','E','F','G','H','I','J','K',
'L','M','N','O','P','Q','R','S','T' ,'U' ,'V' ,'W','Y','Z' };
void ini_DNAconstants(bool gapsRdiffs);

class DNA {
public:
	DNA(string seq, uint SeqL, string names) :isNT(true),Seq(seq), SeqLength(SeqL),
		head(names) {}
	float id(DNA*);
	void rmSeq(uint st, uint en);
	void mask(int st,int en);
	int maskBorderGaps();
	vector<float> createIDprof(vector<DNA*>, options*);
	float fracUnknown();

	bool isNT;
	string Seq;
	uint SeqLength;
	string head;
};

class Ident {
public:
	Ident(vector<DNA*>&, options*);
	void writeAvg(string file);
	void writeDist(string file);
	int maskLowID();
	int rmGapCols();
	void calcID();
	void minGoodPosPerSeq();
	vector<DNA*> getMultiFasta() { return mf; }
private:
	void extend_mask(vector<float>& p, int& st, int& en, float, const string& , int);
	void rmGapsAllSeqs();

	vector<float> avgToOthers;
	vector<vector<float>> disMat;
	//vector<int> gapColCnt;
	//vector<int> gapRowCnt;
	float AvgDist, SdAvgDist;
	bool NT;
	options* opt;
	vector<DNA*> mf;
	vector<uint> stv, env;//temp storage of start / end positions
};
std::istream& safeGetline(std::istream& is, std::string& t);


vector<DNA*> readFasta(string in);
void writeFasta(options*, vector<DNA*>);
