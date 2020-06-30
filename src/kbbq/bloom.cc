#include <stdlib.h>
#include <cstring>
#include <cmath>
#include <random>
#include "bloom.hh"

namespace bloom
{
	std::vector<uint64_t> hash_seq(std::string seq, int k){
		Kmer kmer(k);
		std::vector<uint64_t> ret;
		ret.reserve(seq.length() - k + 1);
		for (int i = 0; i < seq.length(); ++i) {
			kmer.push_back(seq[i]);
			if (kmer.valid()) { //we have a full k-mer
				ret.push_back(kmer.get_hashed());
			}
		}
		return ret;
	}

	void subsample_and_insert(bloomary_t& bfs, std::vector<uint64_t> hashes, double alpha, minion::Random& rng){
		int mask = (1<<PREFIXBITS)-1;
		std::bernoulli_distribution dist(alpha);
		for(uint64_t h: hashes){
			if(dist(rng)){
				bfs[h&mask].insert(h >> PREFIXBITS);
			}
		}
	}

	std::array<std::vector<size_t>,2> overlapping_kmers_in_bf(std::string seq, const bloomary_t& b, int k){
		bloom::Kmer kmer(k);
		std::vector<bool> kmer_present(seq.length()-k+1, false);
		std::vector<size_t> kmers_in(seq.length(), 0);
		std::vector<size_t> kmers_possible(seq.length(), 0);
		size_t incount = 0;
		size_t outcount = 0;
		for(size_t i = 0; i < k; ++i){
			kmer.push_back(seq[i]);
		}
		kmer_present[0] = kmer.valid() ? b[kmer.hashed_prefix()].query(kmer.get_query()) : false;
		for(size_t i = k; i < seq.length(); ++i){
			kmer.push_back(seq[i]);
			kmer_present[i-k+1] = kmer.valid() ? b[kmer.hashed_prefix()].query(kmer.get_query()) : false;
		}
		for(size_t i = 0; i < seq.length(); ++i){
			if(i < seq.length() - k + 1){ //add kmers now in our window
				//debugging
				/* if(seq == "AAGTGGGTTTCTCAGTATTTTATTCTTTTGATATTATCATACATGATACTATCGTCTTGATTTCTTCTTCAGAGAGTTTATTGTTGTTGTAGAAATACAATTGATTTTTGTGTATTGATTTTGTATCCTGCAGCTTTGCTGAATTTTATTT"){
					std::string kmerstr(seq, i, k);
					std::cerr << kmerstr << " " << kmer_present[i] << std::endl;
				} */

				if(kmer_present[i]){
					++incount;
				} else {
					++outcount;
				}
			}
			if(i >= k){ //remove kmers outside our window
				if(kmer_present[i - k]){
					--incount;
				} else {
					--outcount;
				}
			}
			kmers_in[i] = incount;
			kmers_possible[i] = incount + outcount;
		}
		return {kmers_in, kmers_possible};
	}

	int nkmers_in_bf(std::string seq, const bloomary_t& b, int k){
		std::vector<uint64_t> hashes = hash_seq(seq, k);
		int c = 0;
		for(const uint64_t& h: hashes){
			if(b[h&((1<<PREFIXBITS)-1)].query(h>>PREFIXBITS)){
				c++;
			}
		}
		return c;
	}

	std::array<size_t, 2> find_longest_trusted_seq(std::string seq, const bloomary_t& b, int k){
		bloom::Kmer kmer(k);
		size_t anchor_start, anchor_end, anchor_best, anchor_current;
		anchor_start = anchor_end = std::string::npos;
		anchor_best = anchor_current = 0;
		for(size_t i = 0; i < seq.length(); ++i){
			if(kmer.push_back(seq[i]) >= k){
				if(b[kmer.hashed_prefix()].query(kmer.get_query())){
					anchor_current++; //length of current stretch
				}
				else { //we had a streak but the kmer is not trusted
					if(anchor_current > anchor_best){
						anchor_best = anchor_current;
						anchor_end = i - 1; // always > 0 because i >= kmer.size() >= k
						anchor_start = i + 1 - k - anchor_current;
					}
					anchor_current = 0;
				}
			} else if(anchor_current != 0){ //we had a streak but ran into a non-ATGC base
				if(anchor_current > anchor_best){
					anchor_best = anchor_current;
					anchor_end = i - 1; // always > 0 because i >= kmer.size() >= k
					anchor_start = i + 1 - k - anchor_current;
				}
				anchor_current = 0;
			}
		} //we got to the end
		if(anchor_current > anchor_best){
			anchor_best = anchor_current;
			anchor_end = std::string::npos;
			anchor_start = seq.length() + 1 - k - anchor_current;
		}
		return std::array<size_t,2>{{anchor_start, anchor_end}};
	}

	char get_next_trusted_char(const bloom::Kmer& kmer, const bloomary_t& trusted){
		for(char c: {'A','C','G','T'}){
			bloom::Kmer extra = kmer;
			extra.push_back(c);
			if(trusted[extra.hashed_prefix()].query(extra.get_query())){
				return c;
			}
		}
		return 0;
	}

	std::tuple<std::vector<char>,size_t,bool> find_longest_fix(std::string seq, const bloomary_t& trusted, int k){
		// std::cerr << seq << std::endl;
		bloom::Kmer kmer(k);
		std::vector<char> best_c{};
		size_t best_i = 0;
		bool single = false; //this pair of flags will be used to determine whether
		bool multiple = false; //multiple corrections were considered
		char unfixed_char = seq[k-1];
		for(const char& c: {'A','C','G','T'}){
			if(c == unfixed_char){continue;}
			seq[k-1] = c;
			kmer.reset();
			size_t i;
			size_t i_stop = std::max((size_t)2*k-1, seq.length());
			for(i = 0; i < i_stop; ++i){ //i goes to max(2*k-1, seq.length())
				char n = i < seq.length() ? seq[i] : get_next_trusted_char(kmer, trusted);
				if(n == 0){ // no next trusted kmer
					break;
				}
				kmer.push_back(n);
				if(i >= k-1){
					// std::cerr << kmer.size() << std::endl;
					if(kmer.valid()){
#ifndef NDEBUG
						std::cerr << std::string(kmer) << " " << i << " " << trusted[kmer.hashed_prefix()].query(kmer.get_query()) << std::endl;
#endif
						if(!trusted[kmer.hashed_prefix()].query(kmer.get_query())){
							break;
						} else { //we have a trusted kmer with this fix
							if(i == k-1){ //the first possible trusted kmer
								if(single){multiple = true;} //if we had one already, set multiple
								single = true;
							}
						}
					} else { //a non-ATCG base must've been added
						break;
					}
				}
			}
			if(i > best_i){
#ifndef NDEBUG
				std::cerr << std::string(kmer) << " " << i << std::endl;
#endif
				best_c.clear();
				best_c.push_back(c);
				best_i = i;
			} else if (i == best_i){
#ifndef NDEBUG
				std::cerr << std::string(kmer) << " " << i << " Tie" << std::endl;
#endif
				best_c.push_back(c);
			}
		}
		return std::make_tuple(best_c, best_i, multiple);
	}

	Bloom::Bloom(int nshift, int nhashes): ninserts(0),
		//-3 is divide by 8, so we have nshift bits
		bloom(yak_bf_init(std::min(std::max(nshift, YAK_BLK_SHIFT), 64-YAK_BLK_SHIFT), nhashes),
			[](yak_bf_t* p){yak_bf_destroy(p);}) {
		//nshift + YAK_BLK_SHIFT should be less than 64 (nshift <= 55)
		//nshift should be greater than or equal to YAK_BLK_SHIFT (9)
		//thus 9 <= nshift <= 55
	}

	Bloom::~Bloom(){
		// delete this->bloom;
	}

	int Bloom::insert(unsigned long long hash){
		int ninserted = yak_bf_insert(this->bloom.get(), hash);
		if(ninserted < this->bloom->n_hashes){
			this->ninserts++;
		}
		return ninserted;
	}

	int Bloom::query_n(unsigned long long hash) const{
		int x = this->bloom->n_shift - YAK_BLK_SHIFT; // the bloom filter size
		unsigned long long y = hash & ((1ULL<<x)-1); // fit the hash into the bloom filter;
		//discard the beginning (which determines which filter the hash goes into to begin with)
		int h1 = hash >> x & YAK_BLK_MASK; // this is the beginning part that's not in y;
		int h2 = hash >> this->bloom->n_shift & YAK_BLK_MASK; //this is some middle part of the hash;
		uint8_t *p = &this->bloom->b[y << (YAK_BLK_SHIFT-3)];
		if((h2&31) == 0){
			h2 = (h2 + 1) & YAK_BLK_MASK;
		}
		int z = h1;
		int count = 0;
		for(int i = 0; i < this->bloom->n_hashes; z = (z + h2) & YAK_BLK_MASK){
			uint8_t *q = &p[z>>3];
			uint8_t u = 1 << (z&7);
			count += !!(*q & u); // !! = double negate; makes 1 or 0
			++i;
		}
		return count;
	}

	// inline bool Bloom::query(unsigned long long hash) const{
	// 	return (this->query_n(hash) == this->bloom->n_hashes);
	// }

	//fprate given approximate number of times bf was loaded
	double Bloom::fprate(unsigned long long n){
		int m = 1ULL<<(this->bloom->n_shift); //number of bits in the filter
		int k = this->bloom->n_hashes;
		return pow(1 - exp(-1.0 * k * n / m), k);
	}

	//there are 2 ** shift bits in each filter.
	int optimal_nhashes(int shift, int n){
		return ceil(1.0 * (std::pow(2,shift)) / n * log(2));
	}

	long double calculate_fpr(const bloomary_t& bf){
		uint64_t m = 0; // bf[0].nshift; //number of bits
		int k = 0; //bf[0].nhashes; //number of hashes; we take an average and round i thnk
		uint64_t n = 0; //number of insertions
		for(const Bloom& b : bf){
			m += pow(2, b.bloom->n_shift);
			k += b.bloom->n_hashes;
			n += b.ninserts;
		}
		k /= bf.size();
		return pow(1.0l - exp(-(long double)k * (long double)n / (long double)m), (long double)k); // (1 - exp(-kn/m))^k
	}

	long double calculate_phit(const bloomary_t& bf, long double alpha){
		long double fpr = calculate_fpr(bf);
		double exponent = alpha < 0.1 ? 0.2 / alpha : 2;
		long double pa = 1 - pow(1-alpha,exponent);
		return pa + fpr - fpr * pa;
	}

	uint64_t numbits(uint64_t numinserts, long double fpr){
		//m = - n * log2(fpr) / ln2
		return numinserts * (-log2(fpr) / log(2));
	}

	int numhashes(long double fpr){
		//k = -log2(fpr)
		return ceil(-log2(fpr));
	}

	//ensure anchor >= k - 1 before this.
	std::pair<size_t,bool> adjust_right_anchor(size_t anchor, std::string seq, const bloomary_t& trusted, int k){
		bloom::Kmer kmer(k);
		bool multiple = false; //multiple corrections were considered
		size_t modified_idx = anchor + 1;
		assert((anchor >= k - 1));
		for(size_t i = modified_idx - k + 1; i < modified_idx; ++i){
			kmer.push_back(seq[i]);
		}
		for(const char& c : {'A','C','G','T'}){
			if(seq[anchor+1] == c){continue;}
			bloom::Kmer new_kmer = kmer;
			new_kmer.push_back(c);
			//if this fix works, we don't need to adjust the anchor if it fixes all remaining kmers.
			for(size_t i = 0; new_kmer.valid() &&
			trusted[new_kmer.hashed_prefix()].query(new_kmer.get_query()) && modified_idx + 1 + i < seq.length() &&
			i < k; ++i){
#ifndef NDEBUG
				std::cerr << "i: " << i << " anchor + 2 + i: " << anchor + 2 + i << " len: " << seq.length() << std::endl;
#endif
				new_kmer.push_back(seq[modified_idx+1+i]);
				//if we get to the end of the altered kmers and they're all fixed, the anchor is fine.
				if((modified_idx + 1 + i == seq.length()-1 || i == k - 1) && new_kmer.valid() &&
				trusted[new_kmer.hashed_prefix()].query(new_kmer.get_query())){
#ifndef NDEBUG
					std::cerr << "First Try!" << std::endl;
#endif
					return std::make_pair(anchor, multiple);
				}
			}
		} // if we make it through this loop, we need to adjust the anchor.
		//we will test fixes starting with halfway through the last kmer to the end.
		//how much we're winding back the anchor; anchor-i-k must be > 0.
		for(int i = k/2-1; i >= 0 && anchor > i+k; --i){ 
			kmer.reset();
			modified_idx = anchor-i;
			// size_t modified_idx = anchor+1-i+k-2;
			// for(size_t j = anchor + 1 - i - k; j < anchor + 1 - i; ++j){ //start can be -1 from this
			for(size_t j = modified_idx - k + 1; j < modified_idx; ++j){
				kmer.push_back(seq[j]);
			}
			for(const char& c: {'A','C','G','T'}){
				if(seq[modified_idx] == c){continue;}
				bloom::Kmer new_kmer = kmer;
				new_kmer.push_back(c);
#ifndef NDEBUG
				std::cerr << "i: " << i << " modified_idx: " << modified_idx << " kmer:" << seq.substr(modified_idx-k+1,k-1) << c << std::endl;
#endif
				if(new_kmer.valid() && trusted[new_kmer.hashed_prefix()].query(new_kmer.get_query())){
#ifndef NDEBUG
					std::cerr << "Trusted!" << std::endl;
#endif
					multiple = true;
					for(size_t j = 0; new_kmer.valid() &&
					trusted[new_kmer.hashed_prefix()].query(new_kmer.get_query()) &&
					modified_idx+1+j < seq.length() && j <= k/2; ++j){
						new_kmer.push_back(seq[modified_idx+1+j]);
						if(j == k/2 && new_kmer.valid() && trusted[new_kmer.hashed_prefix()].query(new_kmer.get_query())){
							//we went the full length
							return std::make_pair(modified_idx-1, multiple); //the idx before the new base that needs fixing
							//only return here if j == k/2, NOT if you run out of sequence!!!
						}
					}
				}
			}
		}
		//couldn't find a better adjustment
		return std::make_pair(anchor, multiple);
	}

//end namespace
}
