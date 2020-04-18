#pragma once

#include "containers.h"

char DNA_trans[256];
short DNA_amb[256];//to count amb chars
short DNA_IUPAC[256 * 256];
short AA_CODE[256 * 256];
short NT_POS[256];


void ini_DNAconstants(bool gapsRdiffs) {

	short gapPunish = 0;
	if (!gapsRdiffs) { gapPunish = 1; }

	//static char* DNA_trans[256] = "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX";	
		//DNA_trans.resize(256,'X');
	for (int i = 0; i < 256; i++) { DNA_trans[i] = 'N'; }
	DNA_trans['A'] = 'T';	DNA_trans['T'] = 'A';
	DNA_trans['C'] = 'G'; DNA_trans['G'] = 'C';
	DNA_trans['a'] = 'T';	DNA_trans['t'] = 'A';
	DNA_trans['c'] = 'G'; DNA_trans['g'] = 'C';
	for (int i = 0; i < 256; i++) { DNA_amb[i] = 0; }
	DNA_amb['A'] = 1;	DNA_amb['T'] = 1;
	DNA_amb['C'] = 1; DNA_amb['G'] = 1;
	DNA_amb['a'] = 1;	DNA_amb['t'] = 1;
	DNA_amb['c'] = 1; DNA_amb['g'] = 1;


	for (int i = 0; i < 256; i++) { NT_POS[i] = 5; }
	NT_POS['A'] = 0; NT_POS['T'] = 1; NT_POS['G'] = 2; NT_POS['C'] = 3;
	NT_POS['N'] = 4;
	NT_POS['a'] = 0; NT_POS['t'] = 1; NT_POS['g'] = 2; NT_POS['c'] = 3;
	NT_POS['n'] = 4;
	for (int i = 0; i < 256 * 256; i++) { DNA_IUPAC[i] = 1; }
	for (int i = 0; i < 14; i++) {//first: N is always a hit
		DNA_IUPAC['N' + 256 * DNA_SPACE[i]] = 0;
		DNA_IUPAC[256 * 'N' + DNA_SPACE[i]] = 0;
	}
	for (int i = 0; i < 14; i++) {//second: gap is not a hit
		DNA_IUPAC['-' + 256 * DNA_SPACE[i]] = gapPunish;
		DNA_IUPAC[256 * '-' + DNA_SPACE[i]] = gapPunish;
	}
	DNA_IUPAC['-' + 256 * '-'] = 0;
	for (int i = 0; i < 5; i++) {//first: N is always a hit
		DNA_IUPAC['B' + 256 * DNA_SPACE[i]] = 0;
		DNA_IUPAC[256 * 'B' + DNA_SPACE[i]] = 0;
		DNA_IUPAC['H' + 256 * DNA_SPACE[i]] = 0;
		DNA_IUPAC[256 * 'H' + DNA_SPACE[i]] = 0;
		DNA_IUPAC['D' + 256 * DNA_SPACE[i]] = 0;
		DNA_IUPAC[256 * 'D' + DNA_SPACE[i]] = 0;
		DNA_IUPAC['V' + 256 * DNA_SPACE[i]] = 0;
		DNA_IUPAC[256 * 'V' + DNA_SPACE[i]] = 0;
	}

	DNA_IUPAC['B' + 256 * 'A'] = 1; DNA_IUPAC[256 * 'B' + 'A'] = 1;
	DNA_IUPAC[('D' + 256 * 'C')] = 1; DNA_IUPAC[256 * 'D' + 'C'] = 1;
	DNA_IUPAC['H' + 256 * 'G'] = 1; DNA_IUPAC[256 * 'H' + 'G'] = 1;
	DNA_IUPAC['V' + 256 * 'T'] = 1; DNA_IUPAC[256 * 'V' + 'T'] = 1;



	DNA_IUPAC['R' + 256 * 'A'] = 0; DNA_IUPAC[256 * 'R' + 'A'] = 0;
	DNA_IUPAC['R' + 256 * 'G'] = 0; DNA_IUPAC[256 * 'R' + 'G'] = 0;
	DNA_IUPAC['M' + 256 * 'C'] = 0; DNA_IUPAC[256 * 'M' + 'C'] = 0;
	DNA_IUPAC['M' + 256 * 'A'] = 0; DNA_IUPAC[256 * 'M' + 'A'] = 0;
	DNA_IUPAC['Y' + 256 * 'C'] = 0; DNA_IUPAC[256 * 'Y' + 'C'] = 0;
	DNA_IUPAC['Y' + 256 * 'T'] = 0; DNA_IUPAC[256 * 'Y' + 'T'] = 0;
	DNA_IUPAC['K' + 256 * 'G'] = 0; DNA_IUPAC[256 * 'K' + 'G'] = 0;
	DNA_IUPAC['K' + 256 * 'T'] = 0; DNA_IUPAC[256 * 'K' + 'T'] = 0;
	DNA_IUPAC['W' + 256 * 'A'] = 0; DNA_IUPAC[256 * 'W' + 'A'] = 0;
	DNA_IUPAC['W' + 256 * 'T'] = 0; DNA_IUPAC[256 * 'W' + 'T'] = 0;
	DNA_IUPAC['S' + 256 * 'C'] = 0; DNA_IUPAC[256 * 'S' + 'C'] = 0;
	DNA_IUPAC['S' + 256 * 'G'] = 0; DNA_IUPAC[256 * 'S' + 'G'] = 0;

	DNA_IUPAC['A' + 256 * 'A'] = 0; DNA_IUPAC[256 * 'A' + 'A'] = 0;
	DNA_IUPAC['T' + 256 * 'T'] = 0; DNA_IUPAC[256 * 'T' + 'T'] = 0;
	DNA_IUPAC['G' + 256 * 'G'] = 0; DNA_IUPAC[256 * 'G' + 'G'] = 0;
	DNA_IUPAC['C' + 256 * 'C'] = 0; DNA_IUPAC[256 * 'C' + 'C'] = 0;
	//DEBUG
	for (int i = 0; i < 256 * 256; i++) { AA_CODE[i] = 1; }
	for (int i = 0; i < 25; i++) {
		AA_CODE[AA_SPACE[i] + 256 * AA_SPACE[i]] = 0;
		AA_CODE[AA_SPACE[i] + 256 * 'X'] = 0;
		AA_CODE['X' + 256 * AA_SPACE[i]] = 0;
	}
	for (int i = 0; i < 25; i++) {//second: gap is not a hit
		DNA_IUPAC['-' + 256 * AA_SPACE[i]] = gapPunish;
		DNA_IUPAC[256 * '-' + AA_SPACE[i]] = gapPunish;
	}



}

std::istream& safeGetline(std::istream& is, std::string& t)
{
	t.clear();
	//from http://stackoverflow.com/questions/6089231/getting-std-ifstream-to-handle-lf-cr-and-crlf
	// The characters in the stream are read one-by-one using a std::streambuf.
	// That is faster than reading them one-by-one using the std::istream.
	// Code that uses streambuf this way must be guarded by a sentry object.
	// The sentry object performs various tasks,
	// such as thread synchronization and updating the stream state.

	std::istream::sentry se(is, true);
	std::streambuf* sb = is.rdbuf();


	for (;;) {
		int c = sb->sbumpc();
		switch (c) {
		case '\n':
			return is;
		case '\r':
			if (sb->sgetc() == '\n')
				sb->sbumpc();
			return is;
		case EOF:
			// Also handle the case when the last line has no line ending
			if (t.empty())
				is.setstate(std::ios::eofbit);
			return is;
		default:
			t += (char)c;
		}
	}
}
void DNA::rmSeq(uint st, uint en) {
	int diff = en - st;
	if (diff <= 0) { return; }
	Seq.erase(st, diff);
	SeqLength -= diff;
	assert(SeqLength == Seq.length());
}
void DNA::mask(int st, int en) {
	char repl = 'N';
	if (!isNT) { repl = 'X'; }
	if (en > (int)SeqLength) { en = SeqLength; }
	for (int i = st; i < en; i++) {
		Seq[i] = repl;
	}
}
int DNA::maskBorderGaps() {
	
	if (SeqLength <= 0) { return 0; }
	//char repl = 'N';	if (!isNT) { repl = 'X'; }
	int maskedPos = 0;

	bool inGap(Seq[0] == '-'); int st(0);
	bool atStart(true);
	for (uint i = 0; i < SeqLength; i++) {
		if (atStart && inGap ) {
			if (Seq[i] != '-') {
				inGap = false;
				this->mask(st, i);
				maskedPos += i - st;
			}
		}
		else {
			atStart = false;
			if (Seq[i] == '-') {
				if (!inGap) { st = i; inGap = true; }
			}
			else {
				st = 0; inGap = false;
			}
		}
	}
	if (inGap) {
		this->mask(st, SeqLength);
		maskedPos += SeqLength - st;
	}
	return maskedPos;
}


Ident::Ident(vector<DNA*>& mf1, options* opt1):
	avgToOthers(0), disMat(0), AvgDist(0.f), SdAvgDist(0.f),
	opt(opt1),mf(mf1){
	NT = opt->NT;
}
void Ident::calcID() {
	uint numSeq = mf.size();
	bool doSmplAvg = opt->avgToOthers;
	bool doDisM = opt->fullDistMat;
	if (doSmplAvg) {
		avgToOthers.resize(numSeq);
	}
	if (doDisM) {
		disMat = vector < vector <float>>(numSeq, vector<float>(numSeq));
	}
	if (numSeq <= 1) { return; }
	uint cnt(0);
	for (uint i = 0; i < numSeq; i++) {
		//float smplAvg = 0;
		for (uint j = i + 1; j < numSeq; j++) {
			float tmp;
			tmp = mf[i]->id(mf[j]);
			AvgDist += tmp;
			cnt++;
			if (doSmplAvg) {
				avgToOthers[i] += tmp;
				avgToOthers[j] += tmp;

			}
			if (doDisM) {
				disMat[i][j] = tmp;
			}
		}

	}
	if (doSmplAvg) {
		for (uint i = 0; i < numSeq; i++) {
			avgToOthers[i] /= numSeq;
		}
	}
	AvgDist /= cnt;
	for (uint i = 0; i < disMat.size(); ++i) {
		for (uint j = i+1; j < disMat.size(); ++j) {
			SdAvgDist += pow(disMat[i][j] - AvgDist, 2);
		}
	}

	SdAvgDist = sqrt(SdAvgDist / cnt);

}
void Ident::rmGapsAllSeqs() {
	if (stv.size() != env.size()) {
		cerr << "stv != env\n"; exit(6233);
	}
	for (uint j = 0; j < stv.size(); j++) {
		for (uint i = 0; i < mf.size(); i++) {
			mf[i]->rmSeq(stv[j], env[j]);
		}
	}
	//don't forget external call to 
	//stv.clear(); env.clear();
}
int Ident::rmGapCols() {
	int rmGaps(0);
	uint mfs = (uint) mf.size();
	float gfc = opt->gapColFrac;
	if (gfc >= 1.f || mfs==0) { return rmGaps; }
	uint seqL = mf[0]->SeqLength;
	vector<float> gapCnt(seqL,0.f);
	for (uint i = 0; i < mfs; i++) {
		string lS = mf[i]->Seq;
		if (mf[i]->SeqLength != seqL) {
			cerr << "Seq " << mf[i]->head << " seq length " << mf[i]->SeqLength << " != " << seqL<<endl;
			exit(834);
		}
		for (uint j = 0; j < seqL; j++) {
			if (lS[j] == '-') { gapCnt[j] += 1.f; }
		}
	}
	//counted gaps, now calc fraction
	
	int st(0), en(0);
	//clean up vec..
	stv.clear(); env.clear();
	int rmCols(0);
	bool inGap(false);
	for (uint j = 0; j < seqL; j++) {
		gapCnt[j] /= mfs;
		if (gapCnt[j] > gfc) {
			if (!inGap) {
				inGap = true;
				st = j;
			}
		}
		else if (inGap) {
			//rm gap now
			en = j-1;
			inGap = false;
			stv.push_back(st); env.push_back(en);
			rmGaps += en - st;
			en = 0; st = 0;
		}
	}
	rmGapsAllSeqs();
	for (uint i = 0; i < stv.size(); i++) {
		rmCols += env[i] - stv[i];
	}
	stv.clear(); env.clear();
	cerr << "Removed " << rmCols << " columns due to more than " << gfc << " gaps in column\n";

	return rmGaps;
}

void Ident::minGoodPosPerSeq() {
	float gp = opt->minGoodPosFrac;
	if (gp <= 0.f) { return; }
	vector<DNA*> newMF;
	for (uint i = 0; i < mf.size(); i++) {
		if ((1.f-mf[i]->fracUnknown()) < gp) {
			delete mf[i];
		}
		else {
			newMF.push_back(mf[i]);
		}
	}
	int rmN = mf.size() - newMF.size();
	if (rmN > 0) {
		cerr << "Removed "<< rmN<<" sequences due to less than "<<gp<<" good positions\n";
	}
	mf = newMF;
}
void Ident::extend_mask(vector<float>& p, int& st, int& en,float thr,const string& s, int winS) {
	if (en >= (int)p.size()) { en = (int) p.size()-1; }
	int i = (int) en;
	float prevID(p[i]);
	int ps = (uint)p.size();
	int ss = (uint)s.size();
	while (i < ps && (prevID < (p[i] - thr*0.01) || p[i] < thr) ) {
		prevID = p[i]; 
		i++;
	}
	en = i;
	int cnt = 0;
	while (cnt < 10 && i<ps) {
		cnt++;
		if (s[i + winS] == '-') { st = i - 1; break; }
		i++;
		if ((i + winS) >= ss) { break; }
	}


//and extend beginning of gap
	i = (int) st;
	prevID = p[i];
	while (i >=0  && (prevID < (p[i] - thr * 0.01) || p[i] < thr)) {
		prevID = p[i];
		i--;
	}
	st = i;
	//check for gap nearby..
	cnt = 0; 
	while (cnt < 10 && i>=0) {
		cnt++;
		if (s[i + winS] == '-') { st = i+1; break; }
		i--;
	}
}

int Ident::maskLowID() {
	if (!opt->maskLowID) { return 0; }
	if (mf.size() <= 1) { return 0; }
	if (avgToOthers.size() == 0) {
		cerr << "First needs to calc avg dist to other sequences\n";
		exit(3323);
	}
	
	
	//first get expected values
	//float expec1=AvgDist;
	uint maskings(0);
	int winS = opt->winSize;
	int winThr = 50;//how many consecutive NT's at low ID before masking?

	for (uint i = 0; i < mf.size(); i++) {
		//if (avgToOthers[i] > expec1) { continue; }
		bool hasMasked(false);
		vector<float> prof;
		prof = mf[i]->createIDprof(mf, opt);
		//float expec2 = (AvgDist + avgToOthers[i]) / 2;
//		float expec2 = (avgToOthers[i]) - ((1 - AvgDist)*0.25f) - 0.05f;
		float expec2 = (avgToOthers[i]) - ((1 - avgToOthers[i])*0.25f) - max(0.1f,SdAvgDist*3);
		
		//check where a regions falls below this
		int st(-1), en(1); bool inBad(false);
		for (uint y = 0; y < prof.size(); y++) {
			if (inBad) {
				if (prof[y] > expec2) { 
					en = y;
					if (en - st >= winThr) {
						hasMasked = true;
						extend_mask(prof, st, en, avgToOthers[i] - 0.01f,mf[i]->Seq, winS);
						mf[i]->mask(st+ winS, en+ winS);
					}
					//clean up
					st = -1; en = -1;
					inBad = false;
				}
			}
			else {
				if (prof[y] < expec2) {
					st = y;
					inBad = true;
				}
			}
		}
		if (inBad) {
			en = prof.size();
			if (en - st >= winThr) {
				hasMasked = true;
				extend_mask(prof, st, en, avgToOthers[i]-0.01f, mf[i]->Seq,winS);
				mf[i]->mask(st+ winS, en+ winS);
			}
		}
		if (hasMasked) { maskings++; }
	}
	cerr << "Low ID check: Masked " << maskings << " sequences\n";

	return(maskings);
}

void writeFasta(options* opt, vector<DNA*> mf) {
	ofstream out(opt->output);
	if (!out) {
		cerr << "can;t open " << opt->output << endl;
		exit(36);
	}
	for (uint i = 0; i < mf.size(); i++) {
		out << mf[i]->head << endl << mf[i]->Seq << endl;
	}
	out.close();
}

vector<DNA*> readFasta(string inF) {

	istream* in;
	in = new ifstream(inF.c_str());
	if (!(*in)) {
		cerr<<"Cant open file " + inF + "\n";
		exit(11);
	}
	vector<DNA*> ret;
	int cnt(0);
	string tseq(""), head("");
	string line;
	while (safeGetline((*in), line)) {
		//while (getline(in, line, '\n')) {
		cnt++;
		
		if (line[0] != '>') {
			transform(line.begin(), line.end(), line.begin(), ::toupper);
			tseq += line;
		}
		else { // new DNA
			if (tseq.length() > 0) {
				DNA* d = new DNA(tseq, tseq.length(), head);
				if (ret.size() > 0) {
					if (d->SeqLength != ret[0]->SeqLength) {
						cerr << "Seq " << head << " not same seq length as first entry: " << d->SeqLength << " != " << ret[0]->SeqLength << endl;
						exit(02);
					}
				}
				ret.push_back(d);
			}
			head = line;
			tseq = "";
		}
	}
	delete in;
	DNA* d = new DNA(tseq, tseq.length(), head);
	ret.push_back(d);
	return ret;
}

vector<float> DNA::createIDprof(vector<DNA*> mf, options* opt) {
	uint winS = (uint)opt->winSize;
	uint skip = -1;
	vector<float> ret(SeqLength - winS);
	for (uint j = 0; j < mf.size(); j++) {
		if (mf[j]->head == this->head) { skip = j; break; }
	}
	float w_err(0.f), pos(0.f);
	short* code;
	if (isNT) {
		code = DNA_IUPAC;
	}
	else {
		code = AA_CODE;
	}
	for (uint i = 0; i < SeqLength; i++) {
		int i2 = i - winS;
		bool inWindow = i >= winS;
		for (uint j = 0; j < mf.size(); j++) {
			string& oS = mf[j]->Seq;
			if (j == skip) { continue; }
			//don't count unknown sites
			char si = Seq[i]; char osi = oS[i];
			if (si == '-' && osi == '-' ){//|| Seq[i] == 'N' || oS[i] == 'N') {
				if (i > winS) { pos -= 1.f; }
				continue;
			}
			w_err += code[si + 256 * osi];
			if (inWindow) {
				if ((si == '-' && osi == '-')) {//|| Seq[i] == 'N' || oS[i] == 'N') { 
					pos += 1.f; 
				}
				w_err -= code[Seq[i2] + 256 * oS[i2]];
			}
			else {
				pos += 1.f;
			}
		}
		if (i >= winS) {
			ret[i - winS] = ((pos - w_err) / pos);
		}

	}
	return ret;
}

float DNA::fracUnknown() {
	float f = 0.f;
	char nchar = 'N';
	if (!isNT) { nchar = 'X'; }
	for (uint i = 0; i < SeqLength; i++) {
		if (Seq[i] == '-' || Seq[i]==nchar) { f += 1; }
	}
	f /= SeqLength;
	return f;
}

float DNA::id(DNA* o) {
	const string& oS = o->Seq;
	float c_err(0.f), pos(0.f);
	short* code;
	if (isNT) {
		code = DNA_IUPAC;
	}
	else {
		code = AA_CODE;
	}
	for (uint i = 0; i < SeqLength; i++) {
		char si = Seq[i]; char osi = oS[i];
		if (si == '-' && osi == '-') {
			continue;
		}
		//if (Seq[i] != oS[i]) {
		c_err += code[si + 256 * osi];
		//}
		pos += 1.f;
	}
	return((pos - c_err) / pos);
}
