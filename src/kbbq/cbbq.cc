#include "htsiter.hh"
#include "bloom.hh"
#include "covariateutils.hh"
#include "recalibrateutils.hh"
#include <memory>
#include <iostream>
#include <htslib/hfile.h>
#include <getopt.h>

//opens file filename and returns a unique_ptr to the result.
std::unique_ptr<htsiter::HTSFile> open_file(std::string filename, bool is_bam = true){
	std::unique_ptr<htsiter::HTSFile> f(nullptr);
	if(is_bam){
		f = std::move(std::unique_ptr<htsiter::BamFile>(new htsiter::BamFile(filename)));
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
	{"ksize",required_argument,0,'k'},
	{"use-oq",no_argument,0,'u'},
	{"set-oq",no_argument,0,'s'},
	{"genomelen",required_argument,0,'g'},
	{0, 0, 0, 0}
};

int main(int argc, char* argv[]){
	int k = 31;
	uint64_t genomelen = 0; //can est coverage by len / genomelen for fastq
								   //or use index with bam, or estimate given coverage
	uint coverage = 15; //if not given, will be estimated.
	bool set_oq = false;
	bool use_oq = false;

	int opt = 0;
	int opt_idx = 0;
	while((opt = getopt_long(argc,argv,"k:usg:",long_options, &opt_idx)) != -1){
		switch(opt){
			case 'k':
				k = optarg;
				break;
			case 'u':
				use_oq = true;
				break;
			case 's':
				set_oq = true;
				break;
			case 'g':
				genomelen = optarg;
				break;
			case '?':
			default:
				std::cerr << "Unknown argument " << opt << std::endl;
				return 1;
				break;
		}
	}
	if(genomelen <= 0){
		std::cerr << "--genomelen (-g) must be provided and must be > 0" << std::endl;
		return 1;
	}

	std::string filename("-");
	if(optind < argc){
		filename = std::string(argv[optind]);
		while(++optind < argc){
			std::cerr << "Warning: Extra argument " << argv[optind] << " ignored." << std::endl;
		}
	}
	if(coverage = 0){
		//TODO: add ability to specify coverage
		//TODO: check this > 0
	}

	//TODO: flag to use and set OQ flag on BAMs;


	long double alpha = 7.0l / (long double)coverage; // recommended by Lighter authors
	long double sampler_desiredfpr = 0.01;
	long double trusted_desiredfpr = 0.0005;

	//see if we have a bam
	htsFormat fmt;
	hFILE* fp = hopen(filename.c_str(), "r");
	if (hts_detect_format(fp, &fmt) < 0) {
		//error
		std::cerr << "Error opening file " << filename;
		hclose_abruptly(fp);
		return 1;
	}
	hclose(fp);
	bool is_bam = true;
	if(fmt.format == bam || fmt.format == cram){
		is_bam = true;
	} else if (fmt.format == fastq_format){
		is_bam = false;
	} else {
		//error
		std::cerr << "Error: File format must be bam, cram, or fastq.";
		return 1;
	}

	std::unique_ptr<htsiter::HTSFile> file;
	file = std::move(open_file(filename, is_bam));

	htsiter::KmerSubsampler subsampler(file.get(), k, alpha);
	//load subsampled bf x
	recalibrateutils::kmer_cache_t subsampled_hashes = recalibrateutils::subsample_kmers(subsampler);
	bloom::bloomary_t subsampled = init_bloomary(bloom::numbits(genomelen*1.5, sampler_desiredfpr),
		bloom::numhashes(sampler_desiredfpr));
	recalibrateutils::add_kmers_to_bloom(subsampled_hashes, subsampled);

	//calculate thresholds
	//TODO: calculate p any kmer added to bloom filter
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
	file = std::move(open_file(filename, is_bam));
	recalibrateutils::kmer_cache_t trusted_hashes = recalibrateutils::find_trusted_kmers(file.get(), subsampled, thresholds, k);
	bloom::bloomary_t trusted = init_bloomary(bloom::numbits(genomelen*1.5, trusted_desiredfpr),
		bloom::numhashes(trusted_desiredfpr));
	recalibrateutils::add_kmers_to_bloom(trusted_hashes, trusted);

	//use trusted kmers to find errors
	file = std::move(open_file(filename, is_bam));
	covariateutils::CCovariateData data = recalibrateutils::get_covariatedata(file.get(), trusted, k);

	//recalibrate reads and write to file
	covariateutils::dq_t dqs = data.get_dqs();
	file = std::move(open_file(filename, is_bam));
	recalibrateutils::recalibrate_and_write(file.get(), dqs, "-");
	return 0;
}


