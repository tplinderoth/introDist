/*
 * introDist.h
 */

#ifndef INTRODIST_H_
#define INTRODIST_H_

#include <fstream>
#include <utility>
#include <map>

typedef std::pair <unsigned int,double> sitedist; // intended to store position (first) and distance (second)

struct DistCodes {
public:
	DistCodes ();
	std::map<std::string, int> codes;
};

struct WindowInfo {
	unsigned int start;
	unsigned int stop;
	unsigned int nsites;
	double d;
};

void mainInfo ();
void comp2DistInfo (unsigned int windsize, unsigned int stepsize, int useAvg);
int comp2Dist (std::fstream &freqfile, std::fstream &genofile, unsigned int window, unsigned int step, int useAvg, std::string* genofname);
int comp2DistReg (std::fstream &freqfile, std::fstream &genofile, unsigned int window, unsigned int step);
int comp2DistPopAvg (std::fstream &freqfile, std::fstream &genofile, unsigned int window, unsigned int step, std::string* genofname);
int comp2DistIndAvg (std::fstream &freqfile, std::fstream &genofile, unsigned int window, unsigned int step, std::string* genofname);
void comp3DistInfo (unsigned int windsize, unsigned int stepsize);
int parseArgs (int c, char **v, std::fstream &freqfile, std::fstream &genofile, unsigned int &window, unsigned int &step, std::string* genofname, int &useAvg, std::string* dist_type);
int distAnalysis(int argc, char** argv);
double expect2PopDist (double p1, double p2);
double expectIndPopDist (double g, double p);
WindowInfo makeWinInfo (unsigned int start, unsigned int stop, unsigned int nsites, double dist);
int comp2Dist (std::fstream &freqfile, std::fstream &genofile, unsigned int window, unsigned int step, const char* genofname);
int comp3Dist (std::fstream &freqfile, std::fstream &genofile, unsigned int window, unsigned int step);
double calc3Dist (double g1, double g2, double p1, double p2);

#endif /* INTRODIST_H_ */
