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
			if (kmer.size() >= k) { //we have a full k-mer
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
		kmer_present[0] = b[kmer.hashed_prefix()].query(kmer.get_query());
		for(size_t i = k; i < seq.length(); ++i){
			kmer.push_back(seq[i]);
			if(kmer.valid()){
				kmer_present[i-k+1] = b[kmer.hashed_prefix()].query(kmer.get_query());
			} else {
				kmer_present[i-k+1] = false;
			}
		}
		for(size_t i = 0; i < seq.length(); ++i){
			if(i < seq.length() - k + 1){ //add kmers now in our window
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
		for(uint64_t h: hashes){
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

	std::pair<std::vector<char>, int> find_longest_fix(std::string seq, const bloomary_t& t, int k){
		//this function is probably not very efficient.
		bloom::Kmer kmer(k);
		int best_l = 0;
		std::vector<char> best_c;
		char original_c;
		size_t end = std::min(seq.length(), 2UL * k - 1);
		original_c = seq[k-1];
		for(const char* e : {"A","C","G","T"} ){
			char d = *e;
			if( d != original_c ){
				seq[k-1] = d;
				kmer.reset();
				size_t i;
				for(i = 0; i < end; ++i){
					kmer.push_back(seq[i]);
					if (kmer.size() >= k) { // we find a k-mer
						if(!t[kmer.hashed_prefix()].query(kmer.get_hashed()>>PREFIXBITS)){
							break; //if it's not trusted, end here
						}
					} else if (kmer.size() == 0 && i > 0){ //we ran into a non-ATCG character
						break;
					}
				}
				if(i > best_l){
					best_l = i;
					best_c.clear();
					best_c.push_back(d);
				} else if (i == best_l){
					best_c.push_back(d);
				}
			}
		}
		if(best_l == seq.length() && best_c.size() > 1){ //we ran out of kmers to try and we have a tie
			for(const char* e : {"A","C","G","T"} ){
				char d = *e; //the fix character
				if( d != original_c ){
					seq = seq.substr(0, end); //reset str to original length
					size_t original_len = seq.length();
					kmer.reset(); //reset kmer
					seq[k-1] = d;
					size_t i;
					for(i = (end >= k ? end - k : 0); i < 2 * k - 1; ++i){
						int cur_seqlen = seq.length();
						if(i < seq.length()){
							kmer.push_back(seq[i]);
							if(kmer.size() == 0){ //we had a non-N base somehow; this shouldn't happen
								break;
							}
						} else {
							//we have to come up with extra bases.
							for(const char* ex : {"A","C","G","T"}){
								char extra = *ex; //the extension character
								bloom::Kmer extra_kmer = kmer; //copy the previous kmer
								extra_kmer.push_back(extra);
								if(t[extra_kmer.hashed_prefix()].query(extra_kmer.get_query())){
									kmer.push_back(extra); // add character to end
									seq.push_back(extra); //
									break; //try to extend again
								}
							}
						}
						if(cur_seqlen == seq.length()){ //no extra bases were trusted
							break; //exit extension for this fix
						}
					}
					if(i > best_l){ //we extended the read
						best_l = i;
						best_c.clear();
						best_c.push_back(d);
					} else if (i == best_l){
						best_c.push_back(d);
					}
				}
			}
		}
		return std::make_pair(best_c, best_l);
	}

	Bloom::Bloom(int nshift, int nhashes): ninserts(0),
		//-3 is divide by 8, so we have nshift bits
		bloom(yak_bf_init(nshift, nhashes), [](yak_bf_t* p){yak_bf_destroy(p);}) {
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

	inline bool Bloom::query(unsigned long long hash) const{
		return (this->query_n(hash) == this->bloom->n_hashes);
	}

	//fprate given approximate number of times bf was loaded
	double Bloom::fprate(unsigned long long n){
		int m = 1ULL<<(this->bloom->n_shift); //number of bits in the filter
		int k = this->bloom->n_hashes;
		return pow(1 - exp(-1.0 * k * n / m), k);
	}

	//there are 2 ** shift bits in each filter.
	int optimal_nhashes(int shift, int n){
		return floor(1.0 * (std::pow(2,shift)) / n * log(2));
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
		return floor(-log2(fpr));
	}

}






