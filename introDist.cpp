/*
 * introDist.cpp
 * Use genetic distances to look for evidence of introgression
 *
 * TODO:
 * - add option to use genotype probabilities instead of called genotypes
 * - exception handling
 */

#include "introDist.h"
#include "generalUtils.h"
#include <cstring>
#include <iomanip>
#include <cmath>
#include <vector>

int main (int argc, char **argv) {
	int rv = 0;

	if (argc < 2) {
		mainInfo();
	} else {
		rv = distAnalysis(argc, argv);
	}

	return rv;
}

int distAnalysis(int argc, char** argv) {

	int rv = 0;
	std::string dist_t(argv[1]);
	std::fstream freqfile;
	std::fstream genofile;
	const char* genofname = NULL;
	unsigned int window = 1;
	unsigned int step = 1;
	int useAvg = 1; // whether to use genome or window average for comp2Dist

	if ((rv = parseArgs(argc, argv, freqfile, genofile, window, step, genofname, useAvg, &dist_t)) == 0) {
		if (dist_t.compare("comp2Dist") == 0) {
			rv = comp2Dist(freqfile, genofile, window, step, useAvg, genofname);
		}
		else if (dist_t.compare("comp3Dist") == 0) {
			rv = comp3Dist(freqfile, genofile, window, step);
		}
		else {
			rv = -1;
			std::cerr << "\n" << dist_t << " is an invalid function\n";
		}
	} else if (rv == 1) {
		rv = 0;
	}

	if (freqfile.is_open())
		freqfile.close();
	if (genofile.is_open())
		genofile.close();

	return rv;
}

int parseArgs (int c, char **v, std::fstream &freqfile, std::fstream &genofile, unsigned int &window, unsigned int &step, const char* genofname, int &useAvg, std::string* dist_type) {
	int argpos = 2;
	int increment = 0;

	// check distance function
	DistCodes d;

	if (d.codes.find(*dist_type) == d.codes.end()) {
		std::cerr << "\n" << *dist_type << " is an invalid function\n";
		mainInfo();
		return -1;
	}

	switch (d.codes[*dist_type]) {
		case 2 :
			if (c < 6) {
				comp2DistInfo(window, step, useAvg);
				return 1;
			}
			break;
		case 3 :
			if (c < 6) {
				comp3DistInfo(window, step);
				return 1;
			}
			break;
	}


	// parse arguments
	while (argpos < c) {

		if (strcmp(v[argpos], "-freq") == 0) {
				freqfile.open(v[argpos+1], std::ios::in);
				if (! freqfile) {
					std::cerr << "Unable to open allele frequency file " << v[argpos+1] << "\n";
					return -1;
				}
		}

		if (strcmp(v[argpos], "-geno") == 0) {
			genofname = v[argpos+1];
			genofile.open(v[argpos+1], std::ios::in);
			if (! genofile) {
				std::cerr << "Unable to open genotype file " << v[argpos+1] << "\n";
				return -1;
			}
		}

		if (strcmp(v[argpos], "-window") == 0) {
			window = atoi(v[argpos+1]);
			if (window < 0) {
				std::cerr << "Window size must be positive\n";
				return -1;
			}
		}

		if (strcmp(v[argpos], "-step") == 0) {
			step = atoi(v[argpos+1]);
			if (step < 0) {
				std::cerr << "Step size must be positive\n";
				return -1;
			}
		}

		if (strcmp(v[argpos], "-useAverage") == 0) {
			useAvg = atoi(v[argpos+1]);
			switch (useAvg) {
				case 0 :
					break;
				case 1:
					break;
				default :
					std::cerr << "-useAverage can only be 0 (use region) or 1 (use genome-wide average)\n";
					return -1;
			}
		}

		argpos += 2 + increment;
		increment = 0;
	}

	return 0;
}

int comp2Dist (std::fstream &freqfile, std::fstream &genofile, unsigned int window, unsigned int step, int useAvg, const char* genofname) {
	int comp2return = 0;

	switch (useAvg) {
		case 0 :
			std::cerr << "Using window-wide distance between populations\n";
			comp2return = comp2DistReg(freqfile, genofile, window, step);
			break;
		case 1 :
			std::cerr << "Using genome-wide distance between populations\n";
			comp2return = comp2DistAvg(freqfile, genofile, window, step, genofname);
			break;
		default :
			std::cerr << "Need to specify whether to use genome- or window-wide distance between populations in call to comp2Dist\n";
			comp2return = -1;
	}

	return comp2return;
}

int comp2DistAvg (std::fstream &freqfile, std::fstream &genofile, unsigned int window, unsigned int step, const char* genofname) {

	/*
	 * this calculates the difference in the genetic distance between individual i and population2 in a region and
	 * the expected distance between individual i's population and population2 averaged across the genome
	 *
	 * D = [d_R(i, P2) - d_G(Pi, P2)] / max(d_G(Pi, P2),d_R(i, P2))
	 * subscript G denotes average over the genome
	 * subscript R denotes average across region
	 *
	 * at site s:
	 * expected distance between individual i and population 2: E[d(i, P2)] = |g_i - E[g_p1]|/2
	 * expected difference between the distance of individuals i's population and population 2: d(Pi, P2) = |E[g_p1] - E[g_p2]|/2
	 *
	 * d_G(Pi, P2) = 1/N_G * sum_from=1_to=N_G |E[g_p1] - E[g_p2]|/2, where N_G is total number of sites analyzed
	 *
	 * d_R(i, P2) = 1/N_R * sum_from=1_to=N_R |g_i - E[g_p1]|/2, where N_R is the number of sites in the window
	 *
	 * assumes individual i is from population 1
	 *
	 */

	std::vector<sitedist> dist;
	dist.resize(window);
	std::vector<sitedist>::iterator distiter = dist.begin();
	unsigned int lastidx = window-1;
	unsigned int nsites = 0;
	unsigned int dr_nsites = 0;
	double windD = 0.0;
	unsigned int i = 0;
	unsigned int n = window-step;
	double dg = 0.0;
	unsigned int dg_nsites = 0;

	std::stringstream ss;

	std::string(fline);
	std::string(fchr);
	unsigned int fpos;
	double p1f;
	double p2f;

	std::string(gline);
	std::string(gchr);
	unsigned int gpos;
	int geno;

	getline(freqfile, fline);
	getline(genofile, gline);

	// estimate number of lines of input files to set up vector to hold distances
	ss >> gchr;
	unsigned int nlines = filesize(genofname)/(gchr.length()+2*sizeof(int));
	std::vector<WindowInfo> windist;
	windist.reserve(nlines/step + 10);

	// calculate windows
	while (!fline.empty()) {

		// parse frequency file lines
		ss.str(std::string());
		ss.clear();
		ss.str(fline);
		ss >> fchr >> fpos >> p1f >> p2f;

		// parse genotype file lines
		ss.str(std::string());
		ss.clear();
		ss.str(gline);
		ss >> gchr >> gpos >> geno;

		// ensure that files are synched
		if (gchr.compare(fchr) != 0 || gpos != fpos) {
			std::cerr << "Positions in input files differ\n";
			return -1;
		}

		distiter->first = gpos;
		if (geno < 0) {
			// missing data
			distiter->second = -3.0;
		} else {
			// update window information
			distiter->second = expectIndPopDist(geno, p2f);
			// update genome-wide distance
			dg += expect2PopDist(p1f, p2f);
			++dg_nsites;
		}
		++distiter;
		++nsites;

		if (distiter == dist.end()) {

			// calculate window value
			windD = 0.0;
			dr_nsites = 0;
			for (distiter = dist.begin(); distiter != dist.end(); ++distiter) {
				if (distiter->second < -2.0) continue;
				windD += distiter->second;
				++dr_nsites;
			}
			if (dr_nsites > 0) {
				windD /= (double)dr_nsites;
			} else {
				std::cerr << "Warning: window " << gchr << " " << dist[0].first << "," << dist[lastidx].first << " has no sites with data\n";
			}
			windist.push_back(makeWinInfo(dist[0].first, dist[lastidx].first, dr_nsites, windD));

			// prepare vector for new values
			for (i = 0; i < n; ++i) {
				dist[i].first = dist[step+i].first;
				dist[i].second = dist[step+i].second;
			}

			nsites = n;
			distiter = dist.begin() + n;
		}

		getline(freqfile, fline);
		getline(genofile, gline);
	}

	// calculate last window
	if (nsites > n && nsites < window) {
		unsigned int lastwinidx = nsites-1;
		dr_nsites = 0;
		windD = 0.0;
		for (i=0; i < lastwinidx; ++i) {
			if (dist[i].second < -2.0) continue;
			windD += dist[i].second;
			++dr_nsites;
		}
		if (dr_nsites > 0) {
			windD /= (double)dr_nsites;
		} else {
			std::cerr << "Warning: window " << gchr << " " << dist[0].first << "," << dist[lastwinidx].first << " has no sites with data\n";
		}
		windist.push_back(makeWinInfo(dist[0].first, dist[lastwinidx].first, dr_nsites, windD));
	}

	// calculate and print distance statistic
	if (dg_nsites > 0) {
		dg /= (double)dg_nsites;
	} else {
		std::cerr << "There were no sites with complete data, cannot calculate distance statistic\n";
		return 0; // use better exception handling here
	}
	std::cerr << "\nAverage distance between populations: " << dg << "\n\n";

	unsigned int midpoint = 0;
	double Dstat = 0.0;
	double denom;
	for (std::vector<WindowInfo>::iterator winditer = windist.begin(); winditer != windist.end(); ++winditer) {
		denom = dg > winditer->d ? dg : winditer->d;
		if (denom > 0) {
			Dstat = (winditer->d - dg)/denom;
		}
		else { // use exception handling for this
			std::cerr << "Warning: genome and window distances are zero for region " << gchr << " " << winditer->start << "," << winditer->stop << "\n";
			Dstat = 0.0;
		}
		midpoint = (winditer->start + winditer->stop)/2;
		std::cout << gchr << "\t" << winditer->start << "\t" << winditer->stop << "\t" << midpoint << "\t" << Dstat << "\t" << winditer->nsites << "\n";
	}

	return 0;
}

int comp2DistReg (std::fstream &freqfile, std::fstream &genofile, unsigned int window, unsigned int step) {

	/*
	 * this calculates the difference in the genetic distance between individual i and population2
	 * and the expected distance between individual i's population and population2 in a region
	 *
	 * D = [d_R(i, P2) - d_R(Pi, P2)]
	 * subscript R denotes average across region
	 *
	 * at site s:
	 * expected distance between individual i and population 2: E[d(i, P2)] = |g_i - E[g_p1]|/2
	 * expected difference between the distance of individuals i's population and population 2: d(Pi, P2) = |E[g_p1] - E[g_p2]|/2
	 *
	 * d_R(Pi, P2) = 1/N_R * sum_from=1_to=N_R |E[g_p1] - E[g_p2]|/2, where N_R is total number of sites in the window
	 *
	 * d_R(i, P2) = 1/N_R * sum_from=1_to=N_R |g_i - E[g_p1]|/2, where N_R is the number of sites in the window
	 *
	 * assumes individual i is from population 1
	 *
	 */

	std::vector<sitedist> dist;
	dist.resize(window);
	std::vector<sitedist>::iterator distiter = dist.begin();
	unsigned int lastidx = window-1;
	unsigned int nsites = 0;
	unsigned int dr_nsites = 0;
	double windD = 0.0;
	unsigned int midpoint = 0;
	unsigned int i = 0;
	unsigned int n = window-step;

	std::stringstream ss;

	std::string(fline);
	std::string(fchr);
	unsigned int fpos;
	double p1f;
	double p2f;

	std::string(gline);
	std::string(gchr);
	unsigned int gpos;
	int geno;

	getline(freqfile, fline);
	getline(genofile, gline);

	// calculate windows
	while (!fline.empty()) {

		// parse frequency file lines
		ss.str(std::string());
		ss.clear();
		ss.str(fline);
		ss >> fchr >> fpos >> p1f >> p2f;

		// parse genotype file lines
		ss.str(std::string());
		ss.clear();
		ss.str(gline);
		ss >> gchr >> gpos >> geno;

		// ensure that files are synched
		if (gchr.compare(fchr) != 0 || gpos != fpos) {
			std::cerr << "Positions in input files differ\n";
			return -1;
		}

		// update window
		distiter->first=gpos;
		distiter->second = geno < 0 ? -3.0 : expectIndPopDist(geno, p2f) - expect2PopDist(p1f, p2f);
		++distiter;
		++nsites;

		if (distiter == dist.end()) {

			// calculate window value
			windD = 0.0;
			dr_nsites = 0;
			for (distiter = dist.begin(); distiter != dist.end(); ++distiter) {
				if (distiter->second < -2.0) continue;
				windD += distiter->second;
				++dr_nsites;
			}
			if (dr_nsites > 0) {
				windD /= (double)dr_nsites;
			} else {
				std::cerr << "Warning: window " << gchr << " " << dist[0].first << "," << dist[lastidx].first << " has no sites with data\n";
			}
			midpoint = (dist[0].first + dist[lastidx].first)/2;
			std::cout << gchr << "\t" << dist[0].first << "\t" << dist[lastidx].first << "\t" << midpoint << "\t" << windD << "\t" << dr_nsites << "\n";

			// prepare vector for new values
			for (i = 0; i < n; ++i) {
				dist[i].first = dist[step+i].first;
				dist[i].second = dist[step+i].second;
			}

			nsites = n;
			distiter = dist.begin() + n;
		}

		getline(freqfile, fline);
		getline(genofile, gline);
	}

	// print last window
	if (nsites > n && nsites < window) {
		unsigned int lastwinidx = nsites-1;
		dr_nsites = 0;
		windD = 0.0;
		for (i=0; i < lastwinidx; ++i) {
			if (dist[i].second < -2.0) continue;
			windD += dist[i].second;
			++dr_nsites;
		}
		if (dr_nsites > 0) {
			windD /= (double)dr_nsites;
		} else {
			std::cerr << "Warning: window " << gchr << " " << dist[0].first << "," << dist[lastwinidx].first << " has no sites with data\n";
		}
		midpoint = (dist[0].first + dist[lastwinidx].first)/2;
		std::cout << gchr << "\t" << dist[0].first << "\t" << dist[lastwinidx].first << "\t" << midpoint << "\t" << windD << "\t" << dr_nsites << "\n";
	}

	return 0;
}

double expect2PopDist (double p1, double p2) {
	/*
	 * expected distance between two populations
	 */

	// expected genotypes under HWE
	double eg1 = 2.0*p1*(1.0-p1) + 2.0*p1*p1;
	double eg2 = 2.0*p2*(1.0-p2) + 2.0*p2*p2;
	return fabs(eg1 - eg2)/2.0;
}

double expectIndPopDist (double g, double p) {
	/*
	 * expected distance between individual and population
	 */

	// expected genotype of population under HWE
	double eg = 2.0*p*(1.0-p) + 2.0*p*p;
	return fabs(g - eg)/2.0;
}

WindowInfo makeWinInfo (unsigned int start, unsigned int stop, unsigned int nsites, double dist) {
	WindowInfo window = {start, stop, nsites, dist};
	return window;
}

int comp3Dist (std::fstream &freqfile, std::fstream &genofile, unsigned int window, unsigned int step) {

	std::vector<sitedist> dist;
	dist.resize(window);
	std::vector<sitedist>::iterator distiter = dist.begin();
	unsigned int lastidx = window-1;
	unsigned int nsites = 0;
	unsigned int dr_nsites = 0;
	double windD = 0.0;
	double logwinD = 0.0;
	unsigned int midpoint = 0;
	unsigned int i = 0;
	unsigned int n = window-step;

	std::stringstream ss;

	std::string(fline);
	std::string(fchr);
	unsigned int fpos;
	double p1f;
	double p2f;

	std::string(gline);
	std::string(gchr);
	unsigned int gpos;
	int geno1;
	int geno2;

	getline(freqfile, fline);
	getline(genofile, gline);

	// calculate windows
	while (!fline.empty()) {

		// parse frequency file lines
		ss.str(std::string());
		ss.clear();
		ss.str(fline);
		ss >> fchr >> fpos >> p1f >> p2f;

		// parse genotype file lines
		ss.str(std::string());
		ss.clear();
		ss.str(gline);
		ss >> gchr >> gpos >> geno1 >> geno2;

		// ensure that files are synched
		if (gchr.compare(fchr) != 0 || gpos != fpos) {
			std::cerr << "Positions in input files differ\n";
			return -1;
		}

		// update window
		distiter->first = gpos;
		distiter->second = (geno1 < 0 || geno2 < 0) ? -1.0 : calc3Dist(geno1, geno2, p1f, p2f);
		++distiter;
		++nsites;

		if (distiter == dist.end()) {
			windD = 0.0;
			dr_nsites = 0;
			for (distiter = dist.begin(); distiter != dist.end(); ++distiter) {
				if (distiter->second < 0) continue;
				windD += distiter->second;
				++dr_nsites;
			}
			if (dr_nsites > 0) {
				windD /= (double)dr_nsites;
			} else {
				std::cerr << "Warning: window " << gchr << " " << dist[0].first << "," << dist[lastidx].first << " has no sites with data\n";
			}
			midpoint = (dist[0].first + dist[lastidx].first)/2;
			logwinD = -log(1.0-windD);
			if (logwinD == 0.0) logwinD = 0.0; // no negative zero
			std::cout << gchr << "\t" << dist[0].first << "\t" << dist[lastidx].first << "\t" << midpoint << "\t" << logwinD << "\t" << dr_nsites << "\n";

			// prepare vector for new values
			for (i = 0; i < n; ++i) {
				dist[i].first = dist[step+i].first;
				dist[i].second = dist[step+i].second;
			}

			nsites = n;
			distiter = dist.begin() + n;
		}

		getline(freqfile, fline);
		getline(genofile, gline);
	}

	// print last window
	if (nsites > n && nsites < window) {
		unsigned int lastwinidx = nsites-1;
		windD = 0.0;
		dr_nsites = 0;
		for (i=0; i < lastwinidx; ++i) {
			if (dist[i].second < 0) continue;
			windD += dist[i].second;
			++dr_nsites;
		}
		if (dr_nsites > 0) {
			windD /= (double)nsites;
		} else {
			std::cerr << "Warning: window " << gchr << " " << dist[0].first << "," << dist[lastwinidx].first << " has no sites with data\n";
		}
		midpoint = (dist[0].first + dist[lastwinidx].first)/2;
		logwinD = -log(1.0-windD);
		if (logwinD == 0.0) logwinD = 0.0; // no negative zero
		std::cout << gchr << "\t" << dist[0].first << "\t" << dist[lastwinidx].first << "\t" << midpoint << "\t" << logwinD << "\t" << dr_nsites << "\n";
	}

	return 0;
}

double calc3Dist (double g1, double g2, double p1, double p2) {
	/*
	 * this calculates the average, expected distance between (possibly admixed) individuals and their population
	 * weighted by the distance between the distance between the (possibly admixed) individuals
	 *
	 * D_1 = (1-d(i1, i2)) * d({i1,p1}, {i2,p2})
	 *
	 * at site s:
	 * Distance between individuals: d(i1, i2) = |g(i1,s) - g(i2,s)|/2
	 * Average, expected distance between individual and pop: d({i1,p1}, {i1,p2}) = ( E[|g(i1,s) - g(p1,s)|] + E[|g(i2,s) - g(p2,s)|] )/4
	 *
	 * over N sites:
	 * D_N = 1/N * sum_from=1_to=N [ (1-d(i1,i2)) * d({i1,p1}, {i2,p2}) ]
	 */

	// expected genotypes under HWE
	double eg1 = 2.0*p1*(1.0-p1) + 2.0*p1*p1;
	double eg2 = 2.0*p2*(1.0-p2) + 2.0*p2*p2;

	// expected distance between individuals and populations
	double pdist = (fabs(g1-eg1) + fabs(g2-eg2))/4.0;

	// distance between individuals
	double idist = fabs(g1-g2)/2.0;

	// weighted distance statistic
	double d = (1.0-idist)*pdist;

	return (d);
}

DistCodes::DistCodes () {
	codes.insert(std::pair<std::string, int>("comp2Dist", 2));
	codes.insert(std::pair<std::string, int>("comp3Dist", 3));
}

void mainInfo () {
	int w = 12;
	std::cerr << "\nintroDist [function] [-function_inputs]\n"
	<< "\nFunctions:\n"
	<< std::setw(w) << std::left << "comp2Dist" << std::setw(w) << "Find regions where the genetic distance between an individual and a population is less than expected\n"
	<< std::setw(w) << std::left << "comp3Dist" << std::setw(w) << "Distance between 2 individuals and their respective populations, weighted by distance between the individuals\n"
	<< "\n";
}

void comp2DistInfo (unsigned int windsize, unsigned int stepsize, int useAvg) {
	int w = 12;
	std::cerr << "\nintroDist comp2Dist [-freq] [-geno]\n"
	<< "\nInput:\n"
	<< std::setw(w) << std::left << "-freq" << std::setw(w) << "FILE" << "file with (1) chr, (2) position, (3) population1 allele frequency, (4) population2 allele frequency\n"
	<< std::setw(w) << std::left << "-geno" << std::setw(w) << "FILE" << "file with (1) chr, (2) position, (3) individual genotypes\n"
	<< std::setw(w) << std::left << "-window" << std::setw(w) << "INT" << "Window size in number of sites " << "[" << windsize << "]" << "\n"
	<< std::setw(w) << std::left << "-step" << std::setw(w) << "INT" << "step size in number of sites " << "[" << stepsize << "]" << "\n"
	<< std::setw(w) << std::left << "-useAverage" << std::setw(w) << "INT" << "difference between populations is calculated for (0) each window, or (1) averaged across the genome " << "[" << useAvg << "]" << "\n"
	<< "\nIt is assumed that the individual belongs to population 1\n"
	<< "\nOutput:\n"
	<< "(1) chr\n"
	<< "(2) window start\n"
	<< "(3) window stop\n"
	<< "(4) window position midpoint\n"
	<< "(5) genetic distance statistic\n"
	<< "(6) number of sites with data in window\n"
	<< "\n";
}

void comp3DistInfo (unsigned int windsize, unsigned int stepsize) {
	int w = 12;
	std::cerr << "\nintroDist comp3Dist [-freq] [-geno]\n"
	<< "\nInput:\n"
	<< std::setw(w) << std::left << "-freq" << std::setw(w) << "FILE" << "file with (1) chr, (2) position, (3) population1 allele frequency, (4) population2 allele frequency\n"
	<< std::setw(w) << std::left << "-geno" << std::setw(w) << "FILE" << "file with (1) chr, (2) position, (3) individual1 genotypes, (4) individual2 genotypes\n"
	<< std::setw(w) << std::left << "-window" << std::setw(w) << "INT" << "Window size in number of sites " << "[" << windsize << "]" << "\n"
	<< std::setw(w) << std::left << "-step" << std::setw(w) << "INT" << "step size in number of sites " << "[" << stepsize << "]" << "\n"
	<< "\nIt is assumed that individual1 belongs to population1 and that individual2 belongs to population2\n"
	<< "\nOutput:\n"
	<< "(1) chr\n"
	<< "(2) window start\n"
	<< "(3) window stop\n"
	<< "(4) window position midpoint\n"
	<< "(5) genetic distance statistic\n"
	<< "(6) number of sites with data in window\n"
	<< "\n";
}
