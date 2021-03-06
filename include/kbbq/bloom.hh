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
#include <immintrin.h>

#define PREFIXBITS 10
#define KBBQ_MAX_KMER 32

//the number of prefix hashes is 1<<PREFIXBITS - 1

namespace bloom{

class blocked_bloom_filter: public bloom_filter
{
protected:
	typedef unsigned char v32uqi __attribute__ ((__vector_size__ (32))); //activate gcc vectorization
	typedef long long base_type;
	typedef base_type v4di __attribute__ ((__vector_size__ (32)));
	typedef v4di cell_type;
	static constexpr cell_type cell_type_zero = {0ULL,0ULL,0ULL,0ULL}; //used to initialize
	typedef std::unique_ptr<cell_type, std::function<void(cell_type*)>> table_type;
	// typedef std::unique_ptr<unsigned char, std::function<void(unsigned char*)>> table_type;
public:
	static const size_t block_size = 512; //512 bits = 64 bytes
	table_type bit_table_;
	//TODO: ensure table size is a multiple of block_size
	blocked_bloom_filter(): bloom_filter(){}
	blocked_bloom_filter(const bloom_parameters& p){
		projected_element_count_ = p.projected_element_count;
     	inserted_element_count_ = 0;
     	random_seed_ = (p.random_seed * 0xA5A5A5A5) + 1 ;
     	desired_false_positive_probability_ = p.false_positive_probability;
    	salt_count_ = std::max(p.optimal_parameters.number_of_hashes, 2u);
		table_size_ = p.optimal_parameters.table_size;
		//ensure table fits a full block
		table_size_ += (table_size_ % block_size) != 0 ? block_size - (table_size_ % block_size) : 0;
		generate_unique_salt();
		void* ptr = 0;
		int ret = posix_memalign(&ptr, block_size / bits_per_char,
			table_size_ / bits_per_char);
		if(ret != 0){
			throw std::bad_alloc();
		}
		bit_table_ = table_type(static_cast<cell_type*>(ptr),
			[](cell_type* x){free(x);});
		std::uninitialized_fill_n(bit_table_.get(), table_size_ / bits_per_char / sizeof(cell_type),
			cell_type_zero);
	}
	//delete copy ctor
	blocked_bloom_filter(const blocked_bloom_filter&) = delete;
	//delete copy assignment
	blocked_bloom_filter& operator=(const blocked_bloom_filter&) = delete;
	//move ctor
	blocked_bloom_filter(blocked_bloom_filter&& o){
		salt_count_ = std::move(o.salt_count_);
		table_size_ = std::move(o.table_size_);
		bit_table_ = std::move(o.bit_table_);
		salt_ = std::move(o.salt_);
		projected_element_count_ = std::move(o.projected_element_count_);
		inserted_element_count_ = std::move(o.inserted_element_count_);
		random_seed_ = std::move(o.random_seed_);
		desired_false_positive_probability_ = std::move(o.desired_false_positive_probability_);
	}

	//move function
	inline blocked_bloom_filter& operator=(blocked_bloom_filter&& o){
		if(this != &o){
			salt_count_ = std::move(o.salt_count_);
			table_size_ = std::move(o.table_size_);
			bit_table_ = std::move(o.bit_table_);
			salt_ = std::move(o.salt_);
			projected_element_count_ = std::move(o.projected_element_count_);
			inserted_element_count_ = std::move(o.inserted_element_count_);
			random_seed_ = std::move(o.random_seed_);
			desired_false_positive_probability_ = std::move(o.desired_false_positive_probability_);
		}
		return *this;
	}


	inline virtual size_t num_blocks() const{
		assert((table_size_ % block_size) == 0);
		return table_size_ / block_size;
	}

	// inline virtual void compute_indices(const bloom_type& hash, std::size_t& bit_index, std::size_t& bit) const {
	// 	bit_index = hash % block_size; //which bit in the block?
 //      	bit       = bit_index % bits_per_char; // 
	// }

	inline virtual size_t get_block(const bloom_type& hash) const{
		// how lighter does it
		// return (hash % (table_size_ - block_size + 1)) / bits_per_char; 
		// we intead pick an aligned block.
		// hash index * sizeof(block) / sizeof(cell)
		return (hash % num_blocks()) * (block_size / bits_per_char) / (sizeof(cell_type));
		
	}

	//pick the correct vector inside the block and then the correct
	//base type inside the vector.
	inline virtual std::pair<size_t,size_t> get_vector_unit(const size_t& bit_index) const{
		return std::make_pair((bit_index / bits_per_char) / sizeof(cell_type),
			(bit_index / bits_per_char) % (sizeof(cell_type)/sizeof(base_type)));
	}

	inline virtual void insert(const unsigned char* key_begin, const size_t& length){
		size_t bit_index = 0;
		size_t bit = 0;
		//index in table with first byte of block
		size_t block_idx = get_block(hash_ap(key_begin, length, salt_[0]));
		cell_type* block = bit_table_.get() + block_idx;
		size_t vec, unit;
		for(size_t i = 1; i < salt_.size(); ++i){
			compute_indices(hash_ap(key_begin, length, salt_[i]), bit_index, bit);
			//the bit index is out of 512, the bit is out of 8.
			//we advance to the proper vector, then subscript to the proper byte
			std::tie(vec, unit) = get_vector_unit(bit_index);
			(*(block + vec))[unit] |= static_cast<base_type>(1) << bit;
		}
		++inserted_element_count_;
	}

	template <typename T>
	inline void insert(const T& t){
		insert(reinterpret_cast<const unsigned char*>(&t), sizeof(T));
	}

	inline virtual bool contains(const unsigned char* key_begin, const std::size_t length) const {
		size_t bit_index = 0;
		size_t bit = 0;
		//index in table with first byte of block
		size_t block_idx = get_block(hash_ap(key_begin, length, salt_[0])); 
		cell_type* block = bit_table_.get() + block_idx;
		size_t vec, unit;
		for(size_t i = 1; i < salt_.size(); ++i){
			compute_indices(hash_ap(key_begin, length, salt_[i]), bit_index, bit);
			std::tie(vec, unit) = get_vector_unit(bit_index);
			if( ((*(block + vec))[unit] & static_cast<base_type>(1) << bit) !=
				static_cast<base_type>(1) << bit)
			{
				return false;
			}
		}
		return true;
	}

	template <typename T>
	inline bool contains(const T& t) const {
		return contains(reinterpret_cast<const unsigned char*>(&t),static_cast<std::size_t>(sizeof(T)));
	}

	inline double effective_fpp() const {
		long double c = size() / element_count() ;
		long double lambda = block_size / c;
		long double fpp = 0;
		for(int i = 0; i < 3 * lambda; ++i){ //var = lambda, so 3*lambda should include most anything
			long double k = i;
			long double p_block = std::pow(lambda, k) * std::exp(-lambda) / std::tgammal(k+1);
			long double fpr_inner = std::pow(1.0 - std::exp(-1.0 * salt_.size() * i / block_size), 1.0 * salt_.size());
			fpp += p_block * fpr_inner;
		}
		return fpp;
	}

protected:
	inline virtual void compute_indices(const bloom_type& hash, std::size_t& bit_index, std::size_t& bit) const {
		bit_index = hash % block_size; //which bit in the block?
      	bit       = bit_index % (bits_per_char * sizeof(base_type)); //which bit in the base type?
	}
};

class pattern_blocked_bf: public blocked_bloom_filter
{
protected:
	typedef std::unique_ptr<cell_type, std::function<void(cell_type*)>> pattern_type;
	static const size_t num_patterns = 65536; //4MiB
	pattern_type patterns;
public:
	pattern_blocked_bf(): blocked_bloom_filter(){}
	pattern_blocked_bf(const bloom_parameters& p): blocked_bloom_filter(p){
		void* ptr = 0;
		//we have num_patterns patterns, each with size block_size (in bits)
		int ret = posix_memalign(&ptr, block_size / bits_per_char,
			num_patterns * block_size / bits_per_char);
		if(ret != 0){
			throw std::bad_alloc();
		}
		patterns = pattern_type(static_cast<cell_type*>(ptr),
			[](cell_type* x){free(x);});
		std::uninitialized_fill_n(patterns.get(),
			num_patterns * block_size / bits_per_char / sizeof(cell_type),
			cell_type_zero);
		minion::Random rng;
		rng.Seed(random_seed_); //todo: check that seeding is proper for multiple rng instances
		//ie. the rng in kmersubsampler.
		std::uniform_int_distribution<> d(0, block_size-1);

		std::vector<size_t> possible_bits(block_size);
		std::iota(possible_bits.begin(), possible_bits.end(), 0);
		//begin with a fully shuffled array
		std::shuffle(possible_bits.begin(), possible_bits.end(), rng);
		for(size_t i = 0; i < num_patterns; ++i){
			//first index of the block containing the pattern
			size_t block_start_idx = i * (block_size / bits_per_char) / sizeof(cell_type);
			//sample salt_.size() bits
			for(int j = 0; j < salt_.size(); ++j){
				size_t sampled_bit = d(rng,
					std::uniform_int_distribution<>::param_type{j, block_size - 1});
				std::swap(possible_bits[j], possible_bits[sampled_bit]);
			}
			//set the bits in the appropriate pattern
			size_t vec, unit;
			for(size_t j = 0; j < salt_.size(); ++j){
				size_t sampled_bit_number = possible_bits[j]; //the bit out of [0, 511]
				//the index of the correct vector and byte within the block
				std::tie(vec, unit) = get_vector_unit(sampled_bit_number);
				size_t sampled_bit_idx = block_start_idx + vec;
				//set the index of the bit within the char
				(*(patterns.get() + sampled_bit_idx))[unit] |= (static_cast<base_type>(1) << (sampled_bit_number % (sizeof(base_type) * bits_per_char)));
			}
		}
	}
	//delete copy ctor
	pattern_blocked_bf(const pattern_blocked_bf&) = delete;
	//delete copy assignment
	pattern_blocked_bf& operator=(const pattern_blocked_bf&) = delete;
	//move ctor
	pattern_blocked_bf(pattern_blocked_bf&& o):
		blocked_bloom_filter(std::move(o)), patterns(std::move(o.patterns))
		{}

	//move function
	inline pattern_blocked_bf& operator=(pattern_blocked_bf&& o){
		if(this != &o){
			blocked_bloom_filter::operator=(std::move(o));
			patterns = std::move(o.patterns);
		}
		return *this;
	}

	inline virtual size_t get_pattern(const bloom_type& hash) const{
		size_t pattern_number = hash & (num_patterns - 1); //which block
		return pattern_number * (block_size / bits_per_char) / sizeof(cell_type); // how lighter does it
	}

	inline virtual void insert(const unsigned char* key_begin, const size_t& length){
		//index in table with first byte of block
		size_t block = get_block(hash_ap(key_begin, length, salt_[0]));
		size_t pattern = get_pattern(hash_ap(key_begin, length, salt_[1]));
		static_assert(block_size / bits_per_char / sizeof(cell_type) > 0,
			"Block size must be greater than or equal to size of cell type.");
		cell_type* bit_block = reinterpret_cast<cell_type*>(__builtin_assume_aligned(bit_table_.get() + block, block_size / bits_per_char));
		cell_type* pattern_block = reinterpret_cast<cell_type*>(__builtin_assume_aligned(patterns.get() + pattern, block_size / bits_per_char));
		for(size_t i = 0; i < block_size / bits_per_char / sizeof(cell_type); ++i){
			*(bit_block + i) |= *(pattern_block + i);
		}
		++inserted_element_count_;
	}

	template <typename T>
	inline void insert(const T& t)
	{
		// Note: T must be a C++ POD type.
		insert(reinterpret_cast<const unsigned char*>(&t),sizeof(T));
	}

	inline virtual bool contains(const unsigned char* key_begin, const size_t length) const{
		//index in table with first byte of block
		size_t block = get_block(hash_ap(key_begin, length, salt_[0]));
		size_t pattern = get_pattern(hash_ap(key_begin, length, salt_[1]));
		static_assert(block_size / bits_per_char / sizeof(cell_type) > 0,
			"Block size must be greater than or equal to size of cell type.");
		cell_type* bit_block = reinterpret_cast<cell_type*>(__builtin_assume_aligned(bit_table_.get() + block, block_size / bits_per_char));
		cell_type* pattern_block = reinterpret_cast<cell_type*>(__builtin_assume_aligned(patterns.get() + pattern, block_size / bits_per_char));
		for(size_t i = 0; i < block_size / bits_per_char / sizeof(cell_type); ++i){
			if(!_mm256_testc_si256((*(bit_block + i) & *(pattern_block + i)), //==
				*(pattern_block + i)))
			{
				return false;
			}
		}
		return true;
	}

	template <typename T>
	inline bool contains(const T& t) const
	{
		return contains(reinterpret_cast<const unsigned char*>(&t),static_cast<std::size_t>(sizeof(T)));
	}

	inline virtual bloom_type block_hash(const unsigned char* key_begin, const size_t& length) const{
		return hash_ap(key_begin, length, salt_[0]);
	}

	template <typename T>
	inline bloom_type block_hash(const T& t) const{
		return block_hash(reinterpret_cast<const unsigned char*>(&t),static_cast<std::size_t>(sizeof(T)));
	}

	inline virtual bloom_type pattern_hash(const unsigned char* key_begin, const size_t& length) const{
		return hash_ap(key_begin, length, salt_[1]);
	}

	template <typename T>
	inline bloom_type pattern_hash(const T& t) const{
		return pattern_hash(reinterpret_cast<const unsigned char*>(&t),static_cast<std::size_t>(sizeof(T)));
	}

	inline double effective_fpp() const {
		long double c = size() / element_count() ;
		long double lambda = block_size / c;
		long double fpp = 0;
		for(int i = 0; i < 3 * lambda; ++i){ //var = lambda, so 3*lambda should include most anything
			long double p_block = std::pow(lambda, i) * std::exp(-lambda) / std::tgammal(i+1);
			long double p_collision = 1.0l - std::pow(1.0l - 1.0l / (num_patterns), i);
			long double fpr_inner = std::pow(1.0l - std::exp(-1.0l * salt_.size() * i / block_size), 1.0l * salt_.size());
			fpr_inner = p_collision + (1.0l - p_collision) * fpr_inner;
			fpp += p_block * fpr_inner;
		}
		return fpp;
	}

};





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
	typedef pattern_blocked_bf bloom_type;
	Bloom(unsigned long long int projected_element_count, double fpr, unsigned long long int seed = 0xA5A5A5A55A5A5A5AULL);
	Bloom(Bloom&& b) noexcept: bloom(std::move(b.bloom)){} //move ctor
	Bloom& operator=(Bloom&& o){bloom = std::move(o.bloom); return *this;} //move assign
	~Bloom();
	bloom_parameters params;
	bloom_type bloom;
	inline void insert(const Kmer& kmer){if(kmer.valid()){bloom.insert(kmer.get());}}
	template <typename T>
	inline void insert(const T& t){bloom.insert(t);}
	inline bool query(const Kmer& kmer) const {return (kmer.valid() && bloom.contains(kmer.get()));}
	inline double fprate() const {return bloom.effective_fpp();}
	inline unsigned long long inserted_elements() const {return bloom.element_count();}
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