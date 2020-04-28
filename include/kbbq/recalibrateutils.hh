#ifndef RECALIBRATEUTILS_H
#define RECALIBRATEUTILS_H

#include <vector>
#include <string>
#include <cmath>
#include "bloom.hh"
#include "htsiter.hh"
#include "covariateutils.hh"

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

typedef std::array<std::vector<uint64_t>, 1<<PREFIXBITS> kmer_cache_t;

//subsample kmers, hash them, and add them to a prefix tree
//read chunksize kmers at once. if chunksize = -1, read up to MAXINT
kmer_cache_t subsample_kmers(htsiter::KmerSubsampler& s, uint64_t chunksize = -1);

//read all kmers from the subsampler and load them into a cache.
kmer_cache_t read_all_kmers(htsiter::KmerSubsampler& s, uint64_t chunksize = -1);

//add kmers in a given cache to the given assembly of bloom filters.
void add_kmers_to_bloom(kmer_cache_t& kmers, bloom::bloomary_t& filters);

//get some reads from a file, whether a kmer is trusted and put it in a cache.
//this can probably be parallelized if needed because the reads are independent
kmer_cache_t find_trusted_kmers(htsiter::HTSFile* file, bloom::bloomary_t& sampled, std::vector<int> thresholds, int k, uint64_t chunksize = -1);

inline long double q_to_p(int q){return std::pow(10.0l, -((long double)q / 10.0l));}
inline int p_to_q(long double p, int maxscore = 42){return p > 0 ? (int)(-10 * std::log10(p)) : maxscore;}

//get covariate data using the trusted kmers
//this can be parallelized easily since each read is independent
covariateutils::CCovariateData get_covariatedata(htsiter::HTSFile* file, bloom::bloomary_t& trusted, int k);

//recalibrate all reads given the CovariateData
void recalibrate_and_write(htsiter::HTSFile* in, covariateutils::dq_t dqs, std::string outfn);

}

#endif