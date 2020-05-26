#include "htsiter.hh"
#include "bloom.hh"
#include "covariateutils.hh"
#include "recalibrateutils.hh"
#include <memory>
#include <iostream>
#include <htslib/hfile.h>
#include <getopt.h>

//opens file filename and returns a unique_ptr to the result.
std::unique_ptr<htsiter::HTSFile> open_file(std::string filename, bool is_bam = true, bool use_oq = false, bool set_oq = false){
	std::unique_ptr<htsiter::HTSFile> f(nullptr);
	if(is_bam){
		f = std::move(std::unique_ptr<htsiter::BamFile>(new htsiter::BamFile(filename, use_oq, set_oq)));
		// f.reset(new htsiter::BamFile(filename));
	} else {
		f = std::move(std::unique_ptr<htsiter::FastqFile>(new htsiter::FastqFile(filename)));
		// f.reset(new htsiter::FastqFile(filename));
	}
	return f;
}

int check_args(int argc, char* argv[]){
	if(argc < 2){
		std::cerr << "Usage: " << argv[0] << " input.[bam,fq]" << std::endl;
		return 1;
	} else {
		std::cerr << "Selected file: " << std::string(argv[1]) << std::endl;
		return 0;
	}
}

template<typename T>
void print_vec(const std::vector<T>& v){
	std::cerr << "[";
	for(int i = 0; i < v.size()-1; ++i){
		std::cerr << v[i] << ", ";
	}
	std::cerr << v.back() << "]" << std::endl;
}

//nshift is the total number of bits, so it's divided by each block
bloom::bloomary_t init_bloomary(uint64_t nbits, int nhashes){
	bloom::bloomary_t ret;
	int eachbits = nbits / ret.size() + 1; //+1 rounds up
	std::cerr << "Requested a total of " << nbits << " bits." << std::endl;
	std::cerr << "Size of each bloom filter: " << eachbits << " bits." << std::endl;
	std::cerr << "Number of hash functions: " << nhashes << "." << std::endl;
	for(size_t i = 0; i < ret.size(); ++i){
		bloom::Bloom n(log2(eachbits), nhashes);
		ret[i] = std::move(n);
	}
	return ret;
}

//long option, required arg?, flag, value
struct option long_options[] = {
	{"ksize",required_argument,0,'k'}, //default: 31
	{"use-oq",no_argument,0,'u'}, //default: off
	{"set-oq",no_argument,0,'s'}, //default: off
	{"genomelen",required_argument,0,'g'}, //estimated for bam input, required for fastq input
	{"coverage",required_argument,0,'c'}, //default: estimated
	{0, 0, 0, 0}
};

int main(int argc, char* argv[]){
	int k = 31;
	uint64_t genomelen = 0; //est w/ index with bam, w/ fq estimate w/ coverage
	uint coverage = 0; //if not given, will be estimated.
	bool set_oq = false;
	bool use_oq = false;

	int opt = 0;
	int opt_idx = 0;
	while((opt = getopt_long(argc,argv,"k:usg:c:",long_options, &opt_idx)) != -1){
		switch(opt){
			case 'k':
				k = std::stoi(std::string(optarg));
				if(k <= 0 || k > YAK_MAX_KMER){
					std::cerr << "Error: k must be <= " << YAK_MAX_KMER << " and > 0." << std::endl;
				}
				break;
			case 'u':
				use_oq = true;
				break;
			case 's':
				set_oq = true;
				break;
			case 'g':
				genomelen = std::stoull(std::string(optarg));
				break;
			case 'c':
				coverage = std::stoul(std::string(optarg));
				break;
			case '?':
			default:
				std::cerr << "Unknown argument " << opt << std::endl;
				return 1;
				break;
		}
	}

	std::string filename("-");
	if(optind < argc){
		filename = std::string(argv[optind]);
		while(++optind < argc){
			std::cerr << "Warning: Extra argument " << argv[optind] << " ignored." << std::endl;
		}
	}


	long double sampler_desiredfpr = 0.01;
	long double trusted_desiredfpr = 0.0005;

	//see if we have a bam
	htsFormat fmt;
	hFILE* fp = hopen(filename.c_str(), "r");
	if (hts_detect_format(fp, &fmt) < 0) {
		//error
		std::cerr << "Error opening file " << filename << std::endl;
		hclose_abruptly(fp);
		return 1;
	}
	bool is_bam = true;
	if(fmt.format == bam || fmt.format == cram){
		is_bam = true;
	} else if (fmt.format == fastq_format){
		is_bam = false;
	} else {
		//error
		std::cerr << "Error: File format must be bam, cram, or fastq." << std::endl;
		hclose_abruptly(fp);
		return 1;
	}

	if(genomelen == 0){
		if(is_bam){
			std::cerr << "Estimating genome length" << std::endl;
			samFile* sf = hts_hopen(fp, filename.c_str(), "r");
			sam_hdr_t* h = sam_hdr_read(sf);
			for(int i = 0; i < sam_hdr_nref(h); ++i){
				genomelen += sam_hdr_tid2len(h, i);
			}
			sam_hdr_destroy(h);
			hts_close(sf);
			if(genomelen == 0){
				std::cerr << "Header does not contain genome information." <<
					" Unable to estimate genome length; please provide it on the command line" <<
					" using the --genomelen option." << std::endl;
				return 1;
			} else {
				std::cerr << "Genome length is " << genomelen <<" bp." << std::endl;
			}
		} else {
			std::cerr << "Error: --genomelen must be specified if input is not a bam." << std::endl;
		}
	} else {
		if(hclose(fp) != 0){
			std::cerr << "Error closing file!" << std::endl;
		}
	}
	
	std::unique_ptr<htsiter::HTSFile> file;

	if(coverage == 0){
		std::cerr << "Estimating coverage." << std::endl;
		uint64_t seqlen = 0;
		file = std::move(open_file(filename, is_bam, set_oq, use_oq));
		std::string seq("");
		while((seq = file->next_str()) != ""){
			seqlen += seq.length();
		}
		if (seqlen == 0){
			std::cerr << "Error: total sequence length in file " << filename <<
				" is 0. Check that the file isn't empty." << std::endl;
			return 1;
		}
		std::cerr << "Total Sequence length: " << seqlen << std::endl;
		std::cerr << "Genome length: " << genomelen << std::endl;
		coverage = seqlen/genomelen;
		std::cerr << "Estimated coverage: " << coverage << std::endl;
		if(coverage == 0){
			std::cerr << "Error: estimated coverage is 0." << std::endl;
			return 1;
		}
	}

	long double alpha = 7.0l / (long double)coverage; // recommended by Lighter authors
	file = std::move(open_file(filename, is_bam, set_oq, use_oq));

	std::cerr << "Sampling kmers at rate " << alpha << std::endl;

	htsiter::KmerSubsampler subsampler(file.get(), k, alpha);
	//load subsampled bf x
	recalibrateutils::kmer_cache_t subsampled_hashes = recalibrateutils::subsample_kmers(subsampler);

	uint64_t nsampled = 0;
	for(std::vector<uint64_t>& v : subsampled_hashes){
		nsampled += v.size();
	}
	std::cerr << "Sampled " << nsampled << " kmers." << std::endl;

	bloom::bloomary_t subsampled = init_bloomary(bloom::numbits(genomelen*1.5, sampler_desiredfpr),
		bloom::numhashes(sampler_desiredfpr));
	recalibrateutils::add_kmers_to_bloom(subsampled_hashes, subsampled);

	//calculate thresholds
	long double fpr = bloom::calculate_fpr(subsampled);
	uint64_t hits = 0;
	for(bloom::Bloom& b : subsampled){
		hits += b.ninserts;
	}
	std::cerr << "Approximate additions to bloom filter: " << hits << std::endl;
	std::cerr << "Approximate false positive rate: " << fpr << std::endl;
	if(fpr > .15){
		std::cerr << "Error: false positive rate is too high. " <<
			"Increase genomelen parameter and try again." << std::endl;
		return 1;
	}

	long double p = bloom::calculate_phit(subsampled, alpha);
	std::vector<int> thresholds = covariateutils::calculate_thresholds(k, p);

	std::vector<long double> cdf = covariateutils::log_binom_cdf(k,p);

	std::cerr << "Thresholds: [ " ;
	for(auto t : thresholds){std::cerr << t << " ";}
	std::cerr << "]" << std::endl;	

	std::cerr << "log CDF: [ " ;
	for(auto c : cdf){std::cerr << c << " ";}
	std::cerr << "]" << std::endl;

	//get trusted kmers bf using subsampled bf
	std::cerr << "Finding trusted kmers" << std::endl;
	file = std::move(open_file(filename, is_bam, set_oq, use_oq));
	recalibrateutils::kmer_cache_t trusted_hashes = recalibrateutils::find_trusted_kmers(file.get(), subsampled, thresholds, k);
	bloom::bloomary_t trusted = init_bloomary(bloom::numbits(genomelen*1.5, trusted_desiredfpr),
		bloom::numhashes(trusted_desiredfpr));
	recalibrateutils::add_kmers_to_bloom(trusted_hashes, trusted);

	//use trusted kmers to find errors
	std::cerr << "Finding errors" << std::endl;
	file = std::move(open_file(filename, is_bam, set_oq, use_oq));
	covariateutils::CCovariateData data = recalibrateutils::get_covariatedata(file.get(), trusted, k);

	std::vector<std::string> rgvals(readutils::CReadData::rg_to_int.size(), "");
	for(auto i : readutils::CReadData::rg_to_int){
		rgvals[i.second] = i.first;
	}

	std::cerr << "Covariate data:" << std::endl;
	std::cerr << "rgcov:";
	for(int i = 0; i < data.rgcov.size(); ++i){ //rgcov[rg][0] = errors
		std::cerr << i << ": " << rgvals[i] << " {" << data.rgcov[i][0] << ", " << data.rgcov[i][1] << "}" << std::endl;
	}
	std::cerr << "qcov:" << "(" << data.qcov.size() << ")" << std::endl;
	for(int i = 0; i < data.qcov.size(); ++i){
		std::cerr << i << " " << rgvals[i] << "(" << data.qcov[i].size() << ")" << ": [";
		for(int j = 0; j < data.qcov[i].size(); ++j){
			if(data.qcov[i][j][1] != 0){
				std::cerr << j << ":{" << data.qcov[i][j][0] << ", " << data.qcov[i][j][1] << "} ";
			}
		}
		std::cerr << "]" << std::endl;
	}


	//recalibrate reads and write to file
	std::cerr << "Training model" << std::endl;
	covariateutils::dq_t dqs = data.get_dqs();

	std::cerr << "dqs:\n" << "meanq: ";
	print_vec<int>(dqs.meanq);
	std::cerr << "rgdq:";
	print_vec<int>(dqs.rgdq);
	std::cerr << "qscoredq:";
	int i = 0;
	for(const std::vector<int>& v : dqs.qscoredq){
		std::cerr << i++ << ": ";
		print_vec<int>(v);
	}
	// std::cerr << "cycledq:";
	// i = 0;
	// int j = 0;
	// int l = 0;
	// for(const std::vector<std::array<std::vector<int>,2>>& v : dqs.cycledq){
	// 	std::cerr << i++ << ": ";
	// 	for(const std::array<std::vector<int>,2>& a: v){
	// 		std::cerr << j++ << ": ";
	// 		for(const std::vector<int>& b: a){
	// 			std::cerr << l++ << ": ";
	// 			print_vec<int>(b);
	// 		}
	// 	}
	// }
	// i = 0;
	// j = 0;
	// std::cerr << "dinucdq:";
	// for(const std::vector<std::vector<int>> v : dqs.dinucdq){
	// 	std::cerr << i++ << ": ";
	// 	for(const std::vector<int> w : v){
	// 		std::cerr << j++ << ": ";
	// 		print_vec<int>(w);
	// 	}
	// }

	std::cerr << "Recalibrating file" << std::endl;
	file = std::move(open_file(filename, is_bam, set_oq, use_oq));
	recalibrateutils::recalibrate_and_write(file.get(), dqs, "-");
	return 0;
}


