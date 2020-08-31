#ifndef KBBQ_BLOOM_HH
#define KBBQ_BLOOM_HH
#include <cstdint>
#include <utility>
#include <htslib/hts.h>
#include <minion.hpp>
#include <memory>
#include <functional>
#include <iostream>
#include "bloom_filter.hpp"
#include <stdexcept>

#define PREFIXBITS 10
#define KBBQ_MAX_KMER 32

//the number of prefix hashes is 1<<PREFIXBITS - 1

namespace bloom{

//a class to hold an encoded kmer
class Kmer{
protected:
	size_t s; //num times kmer added to since last reset
	int k;
	uint64_t x[2]; //fwd and reverse
	uint64_t mask;
	uint64_t shift;
public:
	Kmer(int k): k(k), s(0), mask(k < 32 ? (1ULL<<k*2) - 1 : -1), shift((k-1)*2) {x[0] = x[1] = 0;}
	Kmer(const Kmer& o): k(o.k), s(o.s), mask(o.mask), shift(o.shift), x{o.x[0], o.x[1]} {}
	//add a character and return the number of times the kmer has been added to since last reset
	inline size_t push_back(char ch){
		int c = seq_nt16_int[seq_nt16_table[ch]];
		if (c < 4){
			x[0] = (x[0] << 2 | c) & mask;                  // forward strand
			x[1] = x[1] >> 2 | (uint64_t)(3 - c) << shift;  // reverse strand
			++s;
		} else this->reset(); // if there is an "N", restart
		return s;
	}
	//get encoded kmer
	inline uint64_t get() const{return x[0] < x[1] ? x[0] : x[1];} //min of x[0] and x[1]
	//get encoded prefix
	inline uint64_t prefix() const{return this->get()&((1<<PREFIXBITS)-1);}
	//empty the kmer and set s to 0
	inline void reset(){s = 0; x[0] = x[1] = 0;}
	//the number of times the kmer has been added to since the last reset
	inline size_t size() const{return s;}
	//whether the kmer has enough bases to be of length k
	inline bool valid() const{return (s >= k);}
	//the length of the kmer
	inline int ksize() const{return k;}
	inline operator std::string() const {
		std::string ret{};
		for(int i = 1; i <= k; ++i){
			ret.push_back( seq_nt16_str[seq_nt16_table['0' + ((x[0] & (3ULL << (2*(k-i)))) >> (2*(k-i)))]] );
		} 
		return ret;
	}
	inline explicit operator bool() const{return this->valid();}
};

//a bloom filter. TODO: make blocked; hold an array of bloom filters and
//delegate insert/query fn's to the appropriate one.
//calculating the right fpr may be a bit involved.
//if we do threads we can also add a mutex for insertion
class Bloom
{
public:
	Bloom(unsigned long long int projected_element_count, double fpr, unsigned long long int seed = 0xA5A5A5A55A5A5A5AULL);
	Bloom(Bloom&& b) noexcept: bloom(std::move(b.bloom)){} //move ctor
	Bloom& operator=(Bloom&& o){bloom = std::move(o.bloom); return *this;} //move assign
	~Bloom();
	bloom_parameters params;
	bloom_filter bloom;
	inline void insert(const Kmer& kmer){if(kmer.valid()){bloom.insert(kmer.get());}}
	template <typename T>
	inline void insert(const T& t){bloom.insert(t);}
	inline bool query(const Kmer& kmer) const {return (kmer.valid() && bloom.contains(kmer.get()));}
	inline double fprate() const {return bloom.effective_fpp();}
	// inline double fprate() const {return bloom.GetActualFP();}
};

// typedef std::array<Bloom,(1<<PREFIXBITS)> bloomary_t;

std::array<std::vector<size_t>,2> overlapping_kmers_in_bf(std::string seq, const Bloom& b, int k = 31);

//return the total number of kmers in b
int nkmers_in_bf(std::string seq, const Bloom& b, int k);

//given a kmer, get the next character (in ACGT order) that would create a trusted
//kmer when appended and return it. Return 0 if none would be trusted.
//Set the flag to test in TGCA order instead.
char get_next_trusted_char(const bloom::Kmer& kmer, const Bloom& trusted, bool reverse_test_order = false);

//return the INCLUSIVE indices bounding the largest stretch of trusted sequence
//if the first value is -1, there are no trusted kmers.
//if the 2nd value is -1 (== std::string::npos), until the end of the string is trusted.
//thus the whole string being trusted looks like {0, std::string::npos}
//while no part of the string being trusted looks like {std::string::npos, std::string::npos}
std::array<size_t,2> find_longest_trusted_seq(std::string seq, const Bloom& b, int k);

//find the longest possible fix for the kmer at position (k-1) until the end
//return the best character (multiple in case of a tie), the index of the next untrusted base,
// and whether multiple corrections were considered for the fix.
//if the length of the fix character vector is 0, no fix was found and correction should end.
//If the sequence is reverse-complemented, set the flag to test bases in reverse order
// (TGCA) instead of (ACGT)
std::tuple<std::vector<char>, size_t, bool> find_longest_fix(std::string seq, const Bloom& t, int k, bool reverse_test_order = false);

//given the sampling rate, calculate the probability any kmer is in the array.
long double calculate_phit(const Bloom& bf, long double alpha);

//given the number of inserts and the desired fpr, calculate the total size of the hash needed
uint64_t numbits(uint64_t numinserts, long double fpr);

//given the desired fpr, calculate the number of hash fn's needed
//assuming we use the optimal number of bits
int numhashes(long double fpr);

//given a sequence and an anchor, see if the anchor could be improved by moving it to the left.
//this is used in the correction step.
//ensure anchor >= k before this function is called.
//return {anchor, multiple}, the location of the new anchor and whether multiple corrections
//were possible during anchor adjustment.
std::pair<size_t, bool> adjust_right_anchor(size_t anchor, std::string seq, const Bloom& trusted, int k);

//get the biggest consecutive trusted block, for creating a trusted anchor when
//one doesn't exist. If there are too many consecutive misses, the procedure
//will end early. 94518
int biggest_consecutive_trusted_block(std::string seq, const Bloom& trusted, int k, int current_len);

}

inline std::ostream& operator<< (std::ostream& stream, const bloom::Kmer& kmer){
	stream << std::string(kmer);
	return stream;
}

#endif