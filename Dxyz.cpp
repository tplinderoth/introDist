/*
 * Dxyz.cpp
 * Author: Tyler Linderoth
 *
 * TODO:
 * 1) Implement better bgzip file reading - Boost only works correctly for gzip
 */

#include "Dxyz.h"
#include <cstring>
#include <iomanip>
#include <cmath>
#include <sstream>
#include <cstdlib>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>

int main (int argc, char** argv) {
	int rv = 0;

	std::ifstream poptable;
	/*
	 * The poptable is a 2-column text file
	 * column 1: sample name in VCF
	 * column 2: INT indicating population, 0 = source, 1 = sink, 2 = outgroup
	*/
	std::ifstream vcf; // input VCF file
	std::ofstream outf; // output file
	double minmaf = 0; // minimum minor allele frequency
	int winsize = 15; // Number of SNPs in a window
	int step = 5; // step size in terms of number of SNPs
	int infmt; // input vcf forma:, 0 = uncompressed, 1 = gzipped/bgzipped, 2 = stdin

	if (!(rv=parseArgs(argc, argv, poptable, vcf, outf, winsize, step, minmaf, infmt))) {
		rv = vcf2Dxyz(poptable, vcf, outf, winsize, step, minmaf, infmt);
	}

	if (poptable.is_open()) poptable.close();
	if (vcf.is_open()) vcf.close();
	if (outf.is_open()) outf.close();

	return rv;
}

void dxyzInfo (int winsize, int step, double maf) {
	int w1 = 12;
	int w2 = 8;

	std::cerr << "\nThis program calculates the Dxyz statistic in sliding windows\n"
	<< "\nTo run:\n"
	<< "Dxzy [arguments]\n\n"
	<< std::setw(w1) << std::left << "-vcf" << std::setw(w2) << std::left << "FILE" << "Input VCF (uncompressed, gzipped/bgzipped), '-' for STDIN\n"
	<< std::setw(w1) << std::left << "-poptable" << std::setw(w2) << std::left << "FILE" << "Population table\n"
	<< std::setw(w1) << std::left << "-out" << std::setw(w2) << std::left << "FILE" << "Output file name [STDOUT]\n"
	<< std::setw(w1) << std::left << "-winsize" << std::setw(w2) << std::left << "INT" << "Number of SNPs in a window [" << winsize << "]\n"
	<< std::setw(w1) << std::left << "-step" << std::setw(w2) << std::left << "INT" << "Window step size in terms of number of SNPs [" << step << "]\n"
	<< std::setw(w1) << std::left << "-maf" << std::setw(w2) << std::left << "FLOAT" << "Minimum MAF among all individuals combined [" << maf << "]\n"
	<< "\nPopulation table format:\n"
	<< "column 1: Sample name from input VCF\n"
	<< "column 2: Value denotating samples populations: 0 = source, 1 = sink, 2 = outgroup\n"
	<< "\nOutput\n"
	<< "(1) Chromosome\n"
	<< "(2) Window start position\n"
	<< "(3) Window end position\n"
	<< "(4) Source-sink Dxy\n"
	<< "(5) Source-outgroup Dxy\n"
	<< "(6) Dxyz statistic\n"
	<< "(7) Average number of sequences with data per SNP\n"
	<< "(8) Number of SNPs used for calculation\n\n";
}

int parseArgs (int argc, char** argv, std::ifstream &poptab, std::ifstream &vcf, std::ofstream &outf, int &winsize, int &step, double &maf, int &infmt) {
	int rv = 0;

		if (argc < 5) {
			dxyzInfo(winsize, step, maf);
			return 1;
		}

		int argpos = 1;
		while (argpos < argc) {
			if (strcmp(argv[argpos], "-vcf") == 0) {
				// input VCF
				const char* invcf_name = argv[argpos+1];
				if (strcmp(invcf_name, "-") == 0) {
					infmt = 2; // reading from standard input
				} else {
					vcf.open(invcf_name, std::ios_base::in | std::ios_base::binary);
					if (! vcf) {
						std::cerr << "Unable to open input VCF " << invcf_name << "\n";
						return -1;
					}
					// check for gzipped file by reading magic numbers
					unsigned char magic [2] = {0};
					vcf.read(reinterpret_cast<char*>(magic), sizeof(magic));
					infmt = (magic[0] == 0x1f && magic[1] == 0x8b) ? 1 : 0;
					vcf.seekg(0, std::ios_base::beg);
				}
			}

			else if (strcmp(argv[argpos], "-poptable") == 0) {
				const char* poptable_name = argv[argpos+1];
				poptab.open(poptable_name, std::ios_base::in);
				if (!poptab) {
					std::cerr << "Unable to open population table " << poptable_name << "\n";
					return -1;
				}
			}

			else if (strcmp(argv[argpos], "-out") == 0) {
				// output file
				const char* outfile = argv[argpos+1];
				outf.open(outfile, std::ios::out);
				if (! outf) {
					std::cerr << "Unable to open output VCF " << *outfile << "\n";
					return -1;
				}
			}

			else if (strcmp(argv[argpos],"-winsize") == 0) {
				// window size
				winsize = atoi(argv[argpos+1]);
				if (winsize < 1) {
					std::cerr << "Window size must be at least one SNP\n";
					return -1;
				}
			}

			else if (strcmp(argv[argpos], "-step") == 0) {
				// step size
				step = atoi(argv[argpos+1]);
				if (step < 1) {
					std::cerr << "Step size must be at least one SNP\n";
					return -1;
				}
			}

			else if (strcmp(argv[argpos], "-maf") == 0) {
				// maf
				maf = atof(argv[argpos+1]);
				if (maf < 0 || maf > 0.5) {
					std::cerr << "-maf must be in range [0,0.5)\n";
					return -1;
				}
			}

			else {
				std::cerr << "Unknown argument " << argv[argpos] << "\n";
				return -1;
			}

			argpos += 2;
		}

		return rv;
}

int vcf2Dxyz (std::ifstream &poptab, std::ifstream &vcf, std::ofstream &outf, int winsize, int step, double minmaf, int infmt) {
		// set up population table map
		std::map <std::string, int> popmap;
		if (setPopTable(poptab, popmap) <= 0) {
			return -1;
		}

		// set up instream
		std::streambuf *inbuf = NULL;
		boost::iostreams::filtering_streambuf<boost::iostreams::input> zipbuf;
		zipbuf.push(boost::iostreams::gzip_decompressor());
		if (infmt == 0) {
			inbuf = vcf.rdbuf();
		}
		else if (infmt == 1) {
			zipbuf.push(vcf);
			inbuf = &zipbuf;
		} else if (infmt == 2) {
			inbuf = std::cin.rdbuf();
		}
		std::istream instream(inbuf);

		// set up outstream
		std::streambuf *outbuf = NULL;
		if (outf.is_open()) {
			outbuf = outf.rdbuf();
		} else {
			outbuf = std::cout.rdbuf();
		}
		std::ostream outstream(outbuf);

		// print header to output
		printHeader(outstream);

		// process VCF file
		std::string vcfline;

		// skip headers and set up vector to hold tokens
		while (getline(instream, vcfline)) {
			if (vcfline[0] != '#' || vcfline[1] != '#')
				break;
		}

		// find individual's positions in allele vector
		std::vector<int> ssidx; // this is the source-sink dxy
		std::vector<int> soidx; // this is the source-outgroup dxy

		std::stringstream ss(vcfline);
		std::string tok;
		int idxpos [2];
		int c = -9;
		int c2 = 0;
		int pop;
		while (ss >> tok) {
			if (c >= 0) {
				if (popmap.find(tok) != popmap.end()) {
					pop = popmap.find(tok)->second;
					idxpos[0] = c2;
					idxpos[1] = c2+1;

					if (pop == 0 || pop == 1) {
							ssidx.insert(ssidx.end(), idxpos, idxpos+2);
							//std::cerr << tok << "\t" << pop << "\t" << ssidx.back()-1 << "\t" << ssidx.back() << "\n"; //debug
					}

					if (pop == 0 || pop == 2) {
							soidx.insert(soidx.end(), idxpos, idxpos+2);
							//std::cerr << tok << "\t" << pop << "\t" << soidx.back()-1 << "\t" << soidx.back() << "\n"; //debug
					}

				} else {
					std::cerr << tok << " is not in population table\n";
					return -1;
				}
				c2 += 2;
			}
			++c;
		}
		unsigned int nind = c;
		int nfields = c+9;
		std::cerr << "Number of individuals in VCF: " << nind << "\n";

		// set up data structures to hold site information

		// matrix to hold Dxy values
		Matrix<double> dxy(winsize, 5); // 0:source-sink dxy, 1:source-out dxy, 2:position, 3:number sequences with data, 4: maf

		// vector to hold alleles
		std::vector<int> alleles;
		alleles.resize(nind*2);

		// loop through VCF sites
		int i = 0, j = 0;
		std::string seqid;
		std::string prevseq;
		int nseqdata; // tracks number of sequences at a site with data
		int nsitedata; // tracks number of sites in a window with data
		double nseqdata_win; // tracks average number of sequences per site in a window with data
		double dxy_components [2]; // 0=Source-Sink Dxy, 1=Source-Outgroup Dxy
		double dxyz;
		double winstart = 0; // holds start position of last printed window
		while (getline(instream, vcfline)) {
		// assume vcf fields are [0]=CHROM, [1]=POS, [2]=ID, [3]=REF, [4]=ALT, [5]=QUAL, [6]=FILTER, [7]=INFO, [8]=FORMAT

			// parse site information
			ss.clear();
			ss.str(vcfline);
			ss >> seqid; // extract CHROM
			ss >> dxy(j,2); // extract POS
			nseqdata = 0;
			i=2;
			while (ss >> tok) {
				if (i > 8) {
					// extract genotype
					if(parseGenotype(tok, alleles, &nseqdata) == alleles.end()) {
						return -1;
					}
				}
				else if (i == 7) {
					dxy(j,4) = getMaf(tok); // extract MAF from INFO
				} else if (i == 8) {
					// ensure genotypes are present
					if (tok[0] != 'G' || tok[1] != 'T') {
						std::cerr << "No GT subfield for " << seqid << " " << dxy(j,2) << "\n";
						return -1;
					}
				}
				++i;
			}

			if (i < nfields) {
				std::cerr << "Input VCF appears truncated\n";
				return -1;
			}

			// calculate average number of pairwise differences for site
			dxy(j,0) = calcDxy(alleles, ssidx); // source-sink dxy
			dxy(j,1) = calcDxy(alleles, soidx); // source-out dxy
			dxy(j,3) = nseqdata; // number of sequences with data
			++j;

			if (j == winsize || (!prevseq.empty() && prevseq != seqid)) {
				winstart=dxy(0,2);

				// calculate Dxyz window
				dxyz = calcDxyzWindow(dxy, dxy.nrow()-1, &nseqdata_win, &nsitedata, minmaf, dxy_components);

				// print window
				printWindow(outstream, dxyz, seqid, dxy(0,2), dxy(winsize-1,2), nseqdata_win, nsitedata, dxy_components);

				// slide dxy window
				j = slideWindow(dxy, winsize, step, j);
			}

			prevseq = seqid;
		}

		// process last window if necessary
		if (dxy(0,2) != winstart) {
			dxyz = calcDxyzWindow(dxy, j-1, &nseqdata_win, &nsitedata, minmaf, dxy_components);
			printWindow(outstream, dxyz, seqid, dxy(0,2), dxy(j-1,2), nseqdata_win, nsitedata, dxy_components);
		}

		return 0;
}

void printHeader (std::ostream &outstream) {
	outstream << "SeqID\tWindowStart\tWindowEnd\tSSDxy\tSODxy\tDxyz\tnSequences\tnSites" << std::endl;
}

void printWindow (std::ostream &outstream, const double &dxyz, const std::string &seqid, const unsigned int &startpos, const unsigned int &endpos,
		const double &nseq, const double &nsites, double* dxy_comp) {
	outstream << seqid << "\t" << startpos << "\t" << endpos << "\t" << dxy_comp[0] << "\t" << dxy_comp[1] << "\t";
	if (dxyz != -9) {
		outstream << dxyz << "\t";
	} else {
		outstream << "NA\t";
	}
	outstream << nseq << "\t" << nsites << "\n";
}

int setPopTable(std::ifstream &popfile, std::map<std::string, int> &popmap) {
	std::string sample;
	int pop;
	while (popfile >> sample) {
		popfile >> pop;
		switch (pop) {
			case 0:
				break;
			case 1:
				break;
			case 2:
				break;
			default:
				std::cerr << "Invalid population label " << pop << "\n";
				return -1;
		}
		popmap.insert(std::pair<std::string, int>(sample, pop));
	}
	if (popmap.size() == 0) {
		std::cerr << "Population map is empty\n";
	}
	return popmap.size();
}

double getMaf (const std::string &info) {
	// returns MAF based on the AF INFO subfield
	static char f [30];
	f[0] = '\0';
	int j=0;
	for (unsigned int i=0; i<info.size(); ++i) {
		if (i+2 < info.size() && info[i]=='A' && info[i+1]=='F' && info[i+2]=='=') {
			i += 3;
			while (i < info.size() && info[i] != ';' && info[i] != ',') {
				f[j] = info[i];
				++j;
				++i;
			}
			f[j] = '\0';
			break;
		}
	}

	if (!f[0]) {
		return -1.0;
	}

	double freq = atof(f);
	return (freq < 0.5 ? freq : 1.0-freq);
}

std::vector<int>::iterator& parseGenotype (const std::string &gt, std::vector<int> &alleles, int* n) {
	// sets the alleles in the allele vector and advances the iterator to the last set position
	static std::vector<int>::iterator it = alleles.begin();

	if (gt[1] == '/' || gt[1] == '|') {
		for (unsigned int i=0; i<3; i+=2) {
			switch (gt[i]) {
				case '0':
					*it = 0;
					++*n;
					break;
				case '1':
					*it = 1;
					++*n;
					break;
				case '.':
					*it = -9;
					break;
				default:
					std::cerr << "Unrecognized allele " << gt[i] << "\n";
					it = alleles.end();
					return it;
			}
			++it;
		}
	} else {
		std::cerr << "Invalid genotype subfield format\n";
		it=alleles.end();
		return it;
	}

	if (it == alleles.end()) {
		it = alleles.begin();
	}

	return it;
}

double calcDxy(const std::vector<int> &alleles, const std::vector<int> &idx) {
	static std::vector<int>::const_iterator it1;
	static std::vector<int>::const_iterator it2;
	double ncompare = 0;
	double ndiff = 0;
	for (it1=idx.begin(); it1!=idx.end()-1; ++it1) {
		if (alleles[*it1] == -9) continue;
		for (it2=it1+1; it2!=idx.end(); ++it2) {
			if (alleles[*it2] == -9) continue;
			if (alleles[*it1] != alleles[*it2]) ++ndiff;
			++ncompare;
		}
	}

	return(ncompare > 0 ? ndiff/ncompare : -9);
}

double calcDxyzWindow(const Matrix<double> &dxy, unsigned int lastidx, double* nseqavg, int *nsites, const double minmaf, double* dxy_comp) {
	/*
	 * calculates Dxyz statistic as average_Dxy(source_pop,sink_pop)/average_Dxy(source_pop,outgroup_pop)
	 * average_Dxy = 1/N*[sum_i=1_to_i=N(Dxy)], N = number of sites
	 *
	 * Encode windows with no data or with no source-outgroup differences as -9
	 */

	dxy_comp[0] = 0; // Source-Sink average dxy
	dxy_comp[1] = 0; // Source-Outgroup average dxy
	*nseqavg = 0;
	*nsites = 0;

	unsigned int i;
	for (i=0; i<=lastidx; ++i) {
		if (dxy(i,0) != -9 && dxy(i,1) != -9 && dxy(i,4) >= minmaf) {
			dxy_comp[0] += dxy(i,0);
			dxy_comp[1] += dxy(i,1);
			*nseqavg += dxy(i,3);
			++*nsites;
		}
	}

	dxy_comp[0] /= *nsites;
	dxy_comp[1] /= *nsites;
	*nseqavg /= *nsites;

	if (*nsites < 1) return -9;

	return (dxy_comp[1] > 0 ? dxy_comp[0]/dxy_comp[1] : -9);
}

unsigned int slideWindow (Matrix<double> &win, const unsigned int winsize, const unsigned int step, unsigned int nsites) {
	/*
	 * Advances window along chromosome
	 * returns the number of sites shifted to the beginning of the window
	 */

	unsigned int nshift = 0;
	if (nsites == winsize) {
		// advance existing window
		nshift = winsize - step;
		for (unsigned int i=0; i<nshift; ++i) {
			int k=step+i;
			for (unsigned int j=0; j<win.ncol(); ++j) {
				win(i,j) = win(k,j);
			}
		}
	}

	return nshift;
}
