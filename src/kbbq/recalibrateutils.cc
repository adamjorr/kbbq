#include "recalibrateutils.hh"

using namespace htsiter;

namespace recalibrateutils{

//if the number of reads doesn't fit in a uint64_t call this twice :)
kmer_cache_t subsample_kmers(KmerSubsampler& s, uint64_t chunksize){
	uint64_t counted = 0;
	kmer_cache_t ret;
	ret.fill(std::vector<uint64_t>());
	//the order here matters since we don't want to advance the iterator if we're chunked out
	for(bloom::Kmer kmer = s.next(); counted++ < chunksize && s.not_eof; kmer = s.next()){
		if(kmer.valid()){
			ret[kmer.prefix()].push_back(kmer.get());
		}
	}
	std::cerr << "Sampled kmers: " << counted << std::endl;
	std::cerr << "Total kmers in dataset: " << s.total_kmers << std::endl;
	return ret;
}

//this is a good target for multithreading ;)
void add_kmers_to_bloom(const kmer_cache_t& kmers, bloom::Bloom& filters){
	for(const std::vector<uint64_t>& v : kmers){
		for(uint64_t kmer : v){
			filters.insert(kmer);
		}
	}
}

kmer_cache_t find_trusted_kmers(HTSFile* file, const bloom::Bloom& sampled, std::vector<int> thresholds, int k, uint64_t chunksize){
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
			kmer.push_back(read.seq[i]);
			if(!read.errors[i]){
				++n_trusted;
			}
			if(i >= k && !read.errors[i-k]){
				--n_trusted;
			}
			if(kmer.valid() && n_trusted == k){
				ret[kmer.prefix()].push_back(kmer.get());
				//trusted kmer here
			}
		}
	}
	return ret;
}

covariateutils::CCovariateData get_covariatedata(HTSFile* file, const bloom::Bloom& trusted, int k){
	covariateutils::CCovariateData data;
#ifndef NDEBUG
	std::ifstream errorsin("../../lighter-debug/corrected.txt");
	int linenum = 0;
	std::string line = "";
#endif
	while(file->next() >= 0){
		readutils::CReadData read = file->get();
		read.get_errors(trusted, k, 6);
#ifndef NDEBUG
		//check that errors are same
		std::getline(errorsin, line);
		linenum++;
		std::vector<bool> lighter_errors(line.length(), false);
		std::transform(line.begin(), line.end(), lighter_errors.begin(),
			[](char c) -> bool {return (c == '1');});
		if( lighter_errors != read.errors){
			std::string message("Line num: " + std::to_string(linenum));
			// std::array<size_t,2> anchors = bloom::find_longest_trusted_seq(read.seq, trusted, k);
			// if(anchors[1] - anchors[0] - k + 1 >= k){ //number of trusted kmers >= k
			// 	anchors[1] = bloom::adjust_right_anchor(anchors[1], read.seq, trusted, k);
			// }
			// std::cerr << "Anchors: [" << anchors[0] << ", " << anchors[1] << "]";
			// std::cerr << " (npos is " << std::string::npos << ")\n";
			std::cerr << message << std::endl << "Errors : " ;
			for(const bool& v : read.errors){
				std::cerr << v;
			}
			std::cerr << std::endl;
			std::cerr << "Lighter: " ;
			for(const bool& v : lighter_errors){
				std::cerr << v;
			}
			std::cerr << std::endl;
			std::cerr << "Seq: " << read.seq << std::endl;
			// bloom::Kmer kmer(k);
			// for(const char& c : std::string("CAGAATAGAAAGATTTATAAATTAAATACTC")){
			// 	std::cerr << c << ":" << seq_nt4_table[c] << ":" << kmer.push_back(c) << "," ;
			// }
			// std::cerr << "Last kmer trusted?" << trusted[kmer.hashed_prefix()].query(kmer.get_query()) << std::endl; 
		}
		assert(lighter_errors == read.errors);
#endif
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