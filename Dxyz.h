/*
 * Dxyz.h
 * Author: Tyler Linderoth
 */

#ifndef DXYZ_H_
#define DXYZ_H_

#include <iostream>
#include <fstream>
#include <assert.h>
#include <map>
#include <vector>

template <class T>
struct Matrix
{
public:
	Matrix (size_t n1, size_t n2)
		:   _nrow(0),
			_ncol(0),
			data(NULL)
	{
		if (n1 > 0 && n2 > 0) {
			_nrow = n1;
			_ncol = n2;
			data = new T* [n1];
			for (unsigned int i = 0; i < n1; ++i) {
				data[i] = new T [n2];
			}
		} else {
			std::cerr << "Unable to initialize matrix with " << n1 << "rows and " << n2 << " columns\n";
		}
	}

	~Matrix () {
		for (unsigned int i = 0; i < _nrow; ++i) {
			delete [] data[i];
		}
		delete [] data;
		data = NULL;
	}

	T& operator() (unsigned int i, unsigned int j) {
		assert(i >= 0 && i < _nrow);
		assert(j >= 0 && j < _ncol);

		return data[i][j];
	}

	T operator() (unsigned int i, unsigned int j) const {
		assert(i >= 0 && i < _nrow);
		assert(j >= 0 && j < _ncol);

		return data[i][j];
	}

	unsigned int nrow() const {
		return _nrow;
	}

	unsigned int ncol() const {
		return _ncol;
	}

private:
	unsigned int _nrow;
	unsigned int _ncol;
	T** data;
};

void dxyzInfo (int winsize, int step, double maf);
int parseArgs (int argc, char** argv, std::ifstream &poptab, std::ifstream &vcf, std::ofstream &outf, int &winsize, int &step, double &maf, int &infmt);
int vcf2Dxyz (std::ifstream &poptab, std::ifstream &vcf, std::ofstream &outf, int winsize, int step, double minmaf, int infmt);
void printHeader (std::ostream &outstream);
void printWindow (std::ostream &outstream, const double &dxyz, const std::string &seqid, const unsigned int &startpos, const unsigned int &endpos,
		const double &nseq, const double &nsites, double* dxy_comp);
int setPopTable(std::ifstream &popfile, std::map<std::string, int> &popmap);
double getMaf (const std::string &info);
std::vector<int>::iterator& parseGenotype (const std::string &gt, std::vector<int> &alleles, int* n);
double calcDxy(const std::vector<int> &alleles, const std::vector<int> &idx);
double calcDxyzWindow(const Matrix<double> &dxy, unsigned int lastidx, double* nseqavg, int *nsites, const double minmaf, double* dxy_comp);
unsigned int slideWindow (Matrix<double> &win, const unsigned int winsize, const unsigned int step, unsigned int nsites);

#endif /* DXYZ_H_ */
