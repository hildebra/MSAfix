//small program to check MSA's for correct alignments, remove empty columns etc
#include "containers.h"
const char* MSAF_ver = "2.0";

int main(int argc, char* argv[])
{

	options* opt = new options(argc, argv);
	ini_DNAconstants(opt->GapsPenalty);


	vector<DNA*> muFa = readFasta(opt->input);
	
	Ident* ID(NULL);
	//create obligatory ID object
	ID = new Ident(muFa, opt);
	//filter columns containing too many gaps (or N's?? - no)
	ID->rmGapCols();
	//Border Gaps in individual Seqs can be removed
	int maskedBorder = 0;
	if (opt->maskBorderGap) {
		for (uint i = 0; i < muFa.size(); i++) {
			int tmpMask = muFa[i]->maskBorderGaps();
			if (tmpMask > 0){maskedBorder++;}
		}
	}
	//at this point calculate identity of seqs
	ID->calcID();
	//mask low id regions in alignments .. likely frame-shifts
	ID->maskLowID();
	//remove sequences that don't have enough good positions left
	ID->minGoodPosPerSeq();

	muFa = ID->getMultiFasta();
	
	writeFasta(opt, muFa);

	delete ID;
	delete opt;
	return (0);
}