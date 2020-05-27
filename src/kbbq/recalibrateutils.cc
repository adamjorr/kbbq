#include "recalibrateutils.hh"

using namespace htsiter;

namespace recalibrateutils{

//if the number of reads doesn't fit in a uint64_t call this twice :)
kmer_cache_t subsample_kmers(KmerSubsampler& s, uint64_t chunksize){
	uint64_t kmer = 0;
	uint64_t counted = 0;
	kmer_cache_t ret;
	ret.fill(std::vector<uint64_t>());
	//the order here matters since we don't want to advance the iterator if we're chunked out
	while(counted++ < chunksize && ((kmer = s.next()) != 0 || !s.readseq.empty())){
		uint64_t prefix = kmer & ((1<<PREFIXBITS)-1);
		ret[prefix].push_back(kmer);
	}
	return ret;
}

//this is a good target for multithreading ;)
void add_kmers_to_bloom(const kmer_cache_t& kmers, bloom::bloomary_t& filters){
	for(int i = 0; i < (1<<PREFIXBITS); ++i){
		for(uint64_t kmer : kmers[i]){
			filters[i].insert(kmer >> PREFIXBITS); //remove the suffix and put into the filter
		}
	}
}

kmer_cache_t find_trusted_kmers(HTSFile* file, const bloom::bloomary_t& sampled, std::vector<int> thresholds, int k, uint64_t chunksize){
	uint64_t counted = 0;
	kmer_cache_t ret;
	ret.fill(std::vector<uint64_t>());
	int n_trusted;
	bloom::Kmer kmer(k);
	//the order here matters since we don't want to advance the iterator if we're chunked out
	while(counted++ < chunksize && file->next() >= 0){
		readutils::CReadData read = file->get();
		read.infer_read_errors(sampled, thresholds, k);
		n_trusted = 0;
		kmer.reset();
		for(int i = 0; i < read.seq.length(); ++i){
			if(!read.errors[i]){
				n_trusted++;
			}
			if(i >= k && !read.errors[i-k]){
				n_trusted--;
			}
			kmer.push_back(read.seq[i]);
			if(kmer.size() >= k && n_trusted == k){
				ret[kmer.hashed_prefix()].push_back(kmer.get_hashed());
			}
		}
	}
	return ret;
}

covariateutils::CCovariateData get_covariatedata(HTSFile* file, const bloom::bloomary_t& trusted, int k){
	covariateutils::CCovariateData data;
	while(file->next() >= 0){
		readutils::CReadData read = file->get();
		read.get_errors(trusted, k);
		data.consume_read(read);
	}
	return data;
}

void recalibrate_and_write(HTSFile* in, const covariateutils::dq_t& dqs, std::string outfn){
	if(in->open_out(outfn) < 0){
		//error!! TODO
		return;
	}
	while(in->next() >= 0){
		readutils::CReadData read = in->get();
		std::vector<int> newquals = read.recalibrate(dqs);
		in->recalibrate(newquals);
		if(in->write() < 0){
			//error! TODO
			return;
		}
	}
}

}