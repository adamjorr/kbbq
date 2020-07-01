#include <stdlib.h>
#include <cstring>
#include <cmath>
#include <random>
#include "bloom.hh"

namespace bloom
{

	Bloom::Bloom(unsigned long long int projected_element_count, double fpr,
	unsigned long long int seed): params(){
		params.projected_element_count = projected_element_count;
		params.false_positive_probability = fpr;
		params.random_seed = seed;
		if(!params){
			throw std::invalid_argument("Error: Invalid bloom filter parameters. \
				Adjust parameters and try again.");
		}
		params.compute_optimal_parameters();
		bloom = bloom_filter(params);
	}

	Bloom::~Bloom(){}

	std::array<std::vector<size_t>,2> overlapping_kmers_in_bf(std::string seq, const Bloom& b, int k){
		bloom::Kmer kmer(k);
		std::vector<bool> kmer_present(seq.length()-k+1, false);
		std::vector<size_t> kmers_in(seq.length(), 0);
		std::vector<size_t> kmers_possible(seq.length(), 0);
		size_t incount = 0;
		size_t outcount = 0;
		for(size_t i = 0; i < k; ++i){
			kmer.push_back(seq[i]);
		}
		kmer_present[0] = kmer.valid() ? b.query(kmer) : false;
		for(size_t i = k; i < seq.length(); ++i){
			kmer.push_back(seq[i]);
			kmer_present[i-k+1] = kmer.valid() ? b.query(kmer) : false;
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

	int nkmers_in_bf(std::string seq, const Bloom& b, int k){
		Kmer kmer(k);
		int count = 0;
		for (size_t i = 0; i < seq.length(); ++i) {
			kmer.push_back(seq[i]);
			if (kmer.valid()) { //we have a full k-mer
				if(b.query(kmer)){
					count++;
				}
			}
		}
		return count;
	}

	char get_next_trusted_char(const Kmer& kmer, const Bloom& trusted){
		for(char c: {'A','C','G','T'}){
			Kmer extra = kmer;
			extra.push_back(c);
			if(trusted.query(extra)){
				return c;
			}
		}
		return 0;
	}

	std::array<size_t, 2> find_longest_trusted_seq(std::string seq, const Bloom& b, int k){
		Kmer kmer(k);
		size_t anchor_start, anchor_end, anchor_best, anchor_current;
		anchor_start = anchor_end = std::string::npos;
		anchor_best = anchor_current = 0;
		for(size_t i = 0; i < seq.length(); ++i){
			if(kmer.push_back(seq[i]) >= k){
				if(b.query(kmer)){
					anchor_current++; //length of current stretch
				} else { //we had a streak but the kmer is not trusted
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

	std::tuple<std::vector<char>,size_t,bool> find_longest_fix(std::string seq, const Bloom& trusted, int k){
		// std::cerr << seq << std::endl;
		Kmer kmer(k);
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
						std::cerr << std::string(kmer) << " " << i << " " << trusted.query(kmer) << std::endl;
#endif
						if(!trusted.query(kmer)){
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
				std::cerr << std::string(kmer) << " L" << std::endl;
#endif
				best_c.clear();
				best_c.push_back(c);
				best_i = i;
			} else if (i == best_i){
#ifndef NDEBUG
				std::cerr << std::string(kmer) << " T" << std::endl;
#endif
				best_c.push_back(c);
			}
		}
		return std::make_tuple(best_c, best_i, multiple);
	}

	long double calculate_phit(const Bloom& bf, long double alpha){
		long double fpr = bf.fprate();
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
	std::pair<size_t,bool> adjust_right_anchor(size_t anchor, std::string seq, const Bloom& trusted, int k){
		Kmer kmer(k);
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
			trusted.query(new_kmer) && modified_idx + 1 + i < seq.length() &&
			i < k; ++i){
#ifndef NDEBUG
				std::cerr << "i: " << i << " anchor + 2 + i: " << anchor + 2 + i << " len: " << seq.length() << std::endl;
#endif
				new_kmer.push_back(seq[modified_idx+1+i]);
				//if we get to the end of the altered kmers and they're all fixed, the anchor is fine.
				if((modified_idx + 1 + i == seq.length()-1 || i == k - 1) && new_kmer.valid() &&
				trusted.query(new_kmer)){
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
				if(new_kmer.valid() && trusted.query(new_kmer)){
#ifndef NDEBUG
					std::cerr << "Trusted!" << std::endl;
#endif
					multiple = true;
					for(size_t j = 0; new_kmer.valid() &&
					trusted.query(new_kmer) &&
					modified_idx+1+j < seq.length() && j <= k/2; ++j){
						new_kmer.push_back(seq[modified_idx+1+j]);
						if(j == k/2 && new_kmer.valid() && trusted.query(new_kmer)){
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
