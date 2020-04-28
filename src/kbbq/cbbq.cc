#include "htsiter.hh"
#include "bloom.hh"
#include "covariateutils.hh"
#include "recalibrateutils.hh"
#include <memory>
#include <iostream>

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
		std::cerr << "Usage: " << argv[0] << " input.bam" << std::endl;
		return 1;
	} else {
		std::cerr << "Selected file: " << std::string(argv[1]) << std::endl;
		return 0;
	}
}

int main(int argc, char* argv[]){
	if (check_args(argc, argv) > 0){
		return 1;
	}
	int k = 31;
	bool is_bam = true;
	std::string filename = std::string(argv[1]);
	//todo: add some logic to see if it's a bam
	std::unique_ptr<htsiter::HTSFile> file;
	file = std::move(open_file(filename, is_bam));

	htsiter::KmerSubsampler subsampler(file.get());
	//load subsampled bf x
	recalibrateutils::kmer_cache_t subsampled_hashes = recalibrateutils::subsample_kmers(subsampler);
	bloom::bloomary_t subsampled{};
	recalibrateutils::add_kmers_to_bloom(subsampled_hashes, subsampled);

	//calculate thresholds
	//TODO: calculate p any kmer added to bloom filter
	double p = 0.05;
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
	bloom::bloomary_t trusted;
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


