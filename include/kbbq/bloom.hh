#ifndef KBBQ_BLOOM_HH
#define KBBQ_BLOOM_HH
#include <cstdint>
#include <utility>
#include <htslib/hts.h>
#include <minion.hpp>
#include <memory>
#include <functional>
#include <iostream>

extern "C"{
	#include <yak.h>
	#include <yak-priv.h>
}

#define PREFIXBITS 10
//the number of prefix hashes is 1<<PREFIXBITS - 1

/* 
	A C++ wrapper around Heng Li's bloom filter in yak:
	https://github.com/lh3/yak
 */
//YAK_BLK_SHIFT= 9
//YAK_BLK_MASK = ((1<<(YAK_BLK_SHIFT)) - 1)= 11111111
//                                       3 =       11
//                                       4 =      100
//                                       7 =      111
//                                       9 =     1001
//                   (1<<(YAK_BLK_SHIFT-3) =   100000
//                                      31 =    11111
// BF_SHIFT is the size of the bloom filter (2 ^ BF_SHIFT bits.)

// TODO: use the approximate number of kmers and a desired false positive rate to
//		 calculate the optimal nshift and nhashes parameter

namespace bloom{

inline uint64_t hash(uint64_t key, uint64_t mask){return yak_hash64(key, mask);}

//a class to hold an encoded kmer
class Kmer{
protected:
	size_t s; //num times kmer added to since last reset
	int k;
	uint64_t x[2]; //fwd and reverse
	uint64_t mask;
	uint64_t shift;
public:
	Kmer(int k): k(k), s(0), mask((1ULL<<k*2) - 1), shift((k-1)*2) {x[0] = x[1] = 0;}
	Kmer(const Kmer& o): k(o.k), s(o.s), mask(o.mask), shift(o.shift), x{o.x[0], o.x[1]} {}
	//add a character and return the number of times the kmer has been added to since last reset
	inline size_t push_back(char ch){
		int c = seq_nt4_table[ch];
		if (c < 4){
			x[0] = (x[0] << 2 | c) & mask;                  // forward strand
			x[1] = x[1] >> 2 | (uint64_t)(3 - c) << shift;  // reverse strand
			++s;
		} else this->reset(); // if there is an "N", restart
		return s;
	}
	//get encoded kmer
	inline uint64_t get(){return x[0] < x[1] ? x[0] : x[1];} //min of x[0] and x[1]
	//get hashed kmer
	inline uint64_t get_hashed(){return hash(this->get(),this->mask);}
	//get encoded prefix
	inline uint64_t prefix(){return this->get()&((1<<PREFIXBITS)-1);}
	//get hashed prefix
	inline uint64_t hashed_prefix(){return this->get_hashed()&((1<<PREFIXBITS)-1);}
	//get hashed kmer without the prefix
	inline uint64_t get_query(){return this->get_hashed() >> PREFIXBITS;}
	//empty the kmer and set s to 0
	inline void reset(){s = 0; x[0] = x[1] = 0;}
	//the number of times the kmer has been added to since the last reset
	inline size_t size(){return s;}
	//whether the kmer has enough bases to be of length k
	inline bool valid(){return (s >= k);}
};

class Bloom
{
public:
	// the constructor approximates yak_bf_init()
	Bloom(): Bloom(0, 4){} //22 = approx 512MB
	// across the 2^10 filters (2^PREFIXBITS), use 2^22 bits each. 2^22 * 2^10 = 2^32 bits = 512 MB
	Bloom(int nshift): Bloom(nshift, 4){}
	Bloom(int nshift, int nhashes);
	Bloom(Bloom&& b) noexcept: bloom(std::move(b.bloom)){} //move ctor
	Bloom& operator=(Bloom&& o){ninserts = o.ninserts; bloom = std::move(o.bloom); return *this;} //move assign
	// this is similar to yak_bf_destroy()
	~Bloom();
	// int nshift; //size of hash; 9 <= nshift <= 55
	// int nhashes; //number of hash functions; 4 seems OK
	uint64_t ninserts;
	std::unique_ptr<yak_bf_t, std::function<void(yak_bf_t*)>> bloom;
	//this is similar to yak_bf_insert()
	int insert(unsigned long long hash); //try to insert the hash. return number of hashes that hit
	int query_n(unsigned long long hash) const; //return the number of hashes that hit but don't insert.
	inline bool query(unsigned long long hash) const {
		return (this->query_n(hash) == this->bloom->n_hashes);
	}
	double fprate(unsigned long long n); //fprate given approximate number of times bf was loaded.
};

int optimal_nhashes(int shift, int n); //calculate the optimal nhashes for a size and number of times loaded.

typedef std::array<Bloom,(1<<PREFIXBITS)> bloomary_t;

//hash each kmer in seq and return a vector of hashes.
//the returned vector is not necessarily of length seq.length() - k + 1
//since non-ATGC characters will not be hashable.
std::vector<uint64_t> hash_seq(std::string seq, int k);

//subsample the hashes and insert to bloom filter b.
// static void subsample_and_insert(Bloom& b, std::vector<uint64_t> hashes, double alpha, minion::Random& rng);

//subsample the hashes and insert into the proper bloom filter based on the prefix.
//htsiter::KmerSubsampler is the prefferred way to do this.
void subsample_and_insert(bloomary_t& bfs, std::vector<uint64_t> hashes, double alpha);

std::array<std::vector<size_t>,2> overlapping_kmers_in_bf(std::string seq, const bloomary_t& b, int k = 31);

//return the total number of kmers in b
int nkmers_in_bf(std::string seq, const bloomary_t& b, int k);

//return the INCLUSIVE indices bounding the largest stretch of trusted sequence
//if the first value is -1, there are no trusted kmers.
//if the 2nd value is -1 (== std::string::npos), until the end of the string is trusted.
//thus the whole string being trusted looks like {0, std::string::npos}
//while no part of the string being trusted looks like {std::string::npos, std::string::npos}
std::array<size_t,2> find_longest_trusted_seq(std::string seq, const bloomary_t& b, int k);

//find the longest possible fix for the kmer at position (k-1) until the end
//return the best character (multiple in case of a tie) and the index of the next untrusted base.
//if the length of the fix character vector is 0, no fix was found and correction should end.
std::pair<std::vector<char>, size_t> find_longest_fix(std::string seq, const bloomary_t& t, int k);

//calculate the false positive rate of the given bloom array.
long double calculate_fpr(const bloomary_t& bf);

//given the sampling rate, calculate the probability any kmer is in the array.
long double calculate_phit(const bloomary_t& bf, long double alpha);

//given the number of inserts and the desired fpr, calculate the total size of the hash needed
uint64_t numbits(uint64_t numinserts, long double fpr);

//given the desired fpr, calculate the number of hash fn's needed
//assuming we use the optimal number of bits
int numhashes(long double fpr);

}
#endif