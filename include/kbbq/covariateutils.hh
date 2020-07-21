#ifndef COVARIATEUTILS_H
#define COVARIATEUTILS_H
#define KBBQ_MAXQ 93

#include <vector>
#include <utility>
#include <cmath>
#include <random>
#include <array>
#include <limits>
#include "readutils.hh"
#include "recalibrateutils.hh"

//fwd declare
namespace readutils{
	class CReadData;
}

namespace covariateutils{

typedef std::vector<int> meanq_t;
typedef std::vector<int> rgdq_t;
typedef std::vector<std::vector<int>> qscoredq_t;
typedef std::vector<std::vector<std::array<std::vector<int>,2>>> cycledq_t;
typedef std::vector<std::vector<std::vector<int>>> dinucdq_t;

struct dq_t
{
	meanq_t meanq;
	rgdq_t rgdq;
	qscoredq_t qscoredq;
	cycledq_t cycledq;
	dinucdq_t dinucdq;
};

typedef std::vector<int> prior1_t; //1 prior for each rg
typedef std::vector<std::vector<int>> prior2_t; //1 prior for each rg -> q pair

class NormalPrior{
public:
	static std::vector<long double> normal_prior;
	NormalPrior(){}
	~NormalPrior(){}
	static long double get_normal_prior(size_t j);
};

/*
def _logpmf(self, x, n, p):
    k = floor(x)
    combiln = (gamln(n+1) - (gamln(k+1) + gamln(n-k+1)))
    return combiln + special.xlogy(k, p) + special.xlog1py(n-k, -p)
*/

inline long double log_binom_pmf(unsigned long long k, unsigned long long n, long double p){
	// k+=1; //+1/+2 correction
	// n+=2;
	long double coefficient = (std::lgamma(n+1) - (std::lgamma(k+1) + std::lgamma(n-k+1)));
	return coefficient + (long double)k * std::log(p) + (long double)(n-k) * std::log1p(-p);
}

inline std::vector<long double> log_binom_cdf(unsigned long long k, long double p){
	std::vector<long double> ret(k+1);
	ret[0] = log_binom_pmf(0,k,p);
	for(int i = 1; i <= k; ++i){
		ret[i] = std::log(std::exp(ret[i-1]) + std::exp(log_binom_pmf(i,k,p)));
	}
	return ret;
}

//return vector of each k (k < n) 
inline std::vector<int> calculate_thresholds(unsigned long long k, long double p, long double quartile = .995l){
	std::vector<int> threshold(k+1,0);
	threshold[0] = 0;
	for(int i = 1; i <= k; ++i){
		std::vector<long double> cdf = log_binom_cdf(i,p);
		for(int j = 0; j < cdf.size(); ++j){
			long double logp = cdf[j];
			if(logp >= std::log(quartile)){
				threshold[i] = j;
				break;
			}
		}
	}
	return threshold;
}

inline bool nt_is_not_n(char c){
	return seq_nt16_int[seq_nt16_table[c]] < 4;
}

//ensure there's no N!!!!
inline int8_t dinuc_to_int(char first, char second){
	int8_t f = seq_nt16_int[seq_nt16_table[first]];
	int8_t s = seq_nt16_int[seq_nt16_table[second]];
	return 15 & ((f << 2) | s); // 1111 & (xx00 | 00xx)
}

inline std::array<char, 2> int_to_dinuc(int8_t dinuc){
	std::array<char,2> x = {std::to_string(dinuc >> 2)[0], std::to_string(dinuc | 3)[0]};
	return x;
}

typedef std::array<unsigned long long, 2> covariate_t;

class CCovariate: public std::vector<covariate_t>
{
public:
	CCovariate(): std::vector<covariate_t>(){}
	CCovariate(size_t len): std::vector<covariate_t>(len) {}
	void increment(size_t idx, covariate_t value);
	void increment(size_t idx, unsigned long long err, unsigned long long total);
};

class CRGCovariate: public CCovariate
{
public:
	CRGCovariate(){}
	CRGCovariate(size_t len): CCovariate(len){}
	void consume_read(const readutils::CReadData& read);
	rgdq_t delta_q(prior1_t prior);
};

class CQCovariate: public std::vector<CCovariate>
{
public:
	CQCovariate(){}
	CQCovariate(size_t rgs, size_t qlen): std::vector<CCovariate>(rgs, CCovariate(qlen)){}
	void consume_read(const readutils::CReadData& read);
	qscoredq_t delta_q(prior1_t prior);
};

typedef std::array<CCovariate,2> cycle_t;

//The first Covariate is for fwd reads, the 2nd Covariate is for reverse reads.
class CCycleCovariate: public std::vector<std::vector<cycle_t>>
{
public:
	CCycleCovariate(){}
	CCycleCovariate(size_t rgs, size_t qlen, size_t cylen):
		std::vector<std::vector<cycle_t>>(rgs, std::vector<cycle_t>(qlen, cycle_t({CCovariate(cylen),CCovariate(cylen)})))
		{}
	void consume_read(const readutils::CReadData& read);
	cycledq_t delta_q(prior2_t prior);
};

class CDinucCovariate: public std::vector<std::vector<CCovariate>>
{
public:
	CDinucCovariate(){}
	CDinucCovariate(size_t rgs, size_t qlen, size_t dilen):
		std::vector<std::vector<CCovariate>>(rgs, std::vector<CCovariate>(qlen, CCovariate(dilen)))
		{}
	void consume_read(const readutils::CReadData& read, int minscore = 6);
	dinucdq_t delta_q(prior2_t prior);
};

class CCovariateData
{
// protected:
// 	CRGCovariate rgcov;
// 	CQCovariate qcov;
// 	CCycleCovariate cycov;
// 	CDinucCovariate dicov;
public:
	CRGCovariate rgcov;
	CQCovariate qcov;
	CCycleCovariate cycov;
	CDinucCovariate dicov;
	CCovariateData(){};
	void consume_read(readutils::CReadData& read, int minscore = 6);
	dq_t get_dqs();
};

}
#endif