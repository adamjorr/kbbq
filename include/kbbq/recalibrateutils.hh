#ifndef RECALIBRATEUTILS_H
#define RECALIBRATEUTILS_H

#include <vector>
#include <string>
#include <cmath>
#include "bloom.hh"
#include "htsiter.hh"
#include "covariateutils.hh"

//debug
#include <fstream>


//fwd declare covariateutils stuff
namespace covariateutils{
	class CCovariateData;
	struct dq_t;
}

namespace htsiter{
	class KmerSubsampler;
	class HTSFile;
}

namespace recalibrateutils{

//subsample kmers, hash them, and add them to the bloom filter
void subsample_kmers(htsiter::KmerSubsampler& s, bloom::Bloom& sampled);

//get some reads from a file, whether a kmer is trusted and put it in a cache.
//this can probably be parallelized if needed because the reads are independent
void find_trusted_kmers(htsiter::HTSFile* file, bloom::Bloom& trusted,
	const bloom::Bloom& sampled, std::vector<int> thresholds, int k);

inline long double q_to_p(int q){return std::pow(10.0l, -((long double)q / 10.0l));}
inline int p_to_q(long double p, int maxscore = 42){return p > 0 ? (int)(-10 * std::log10(p)) : maxscore;}

//get covariate data using the trusted kmers
//this can be parallelized easily since each read is independent
covariateutils::CCovariateData get_covariatedata(htsiter::HTSFile* file, const bloom::Bloom& trusted, int k);

//recalibrate all reads given the CovariateData
void recalibrate_and_write(htsiter::HTSFile* in, const covariateutils::dq_t& dqs, std::string outfn);

}

#endif