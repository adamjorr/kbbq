#include <bitset>

int test(){
	bloom_parameters p{};
	p.projected_element_count = 100;
	p.false_positive_probability = .05;
	p.compute_optimal_parameters();
	bloom::pattern_blocked_bf bf(p);
	std::cerr << "Table size: " << bf.size() << std::endl;
	std::cerr << "Block size: " << bf.block_size << std::endl;
	std::cerr << "Block bytes: " << bf.block_size / bits_per_char << std::endl;
	std::cerr << "Pattern bytes: " << bf.num_patterns * bf.block_size / bits_per_char << std::endl;
	std::cerr << "Num hashes: " << bf.hash_count() << std::endl;
	std::cerr << "Optimal K: " << p.optimal_parameters.number_of_hashes << std::endl;
	for(size_t i = 0; i < 10; ++i){
		// std::cerr << (void*)bf.patterns.get() << std::endl;
		std::cerr << i << "(" << (void*)&bf.patterns.get()[2 * i] << ")" << ": ";
		for(size_t j = 0; j < 4; ++j){
			std::cerr << std::bitset<64>(bf.patterns.get()[2*i][j]) << ", ";
		}
		for(size_t j = 0; j < 4; ++j){
			std::cerr << std::bitset<64>(bf.patterns.get()[2*i+1][j]) << ", ";
		}
		std::cerr << std::endl;
	}
	uint64_t kmer = 21345423534512;
	uint64_t gibberish = 123543451243ull;
	bf.insert(kmer);
	std::cerr << "kmer inserted into block: " << bf.get_block(bf.block_hash(kmer)) << std::endl;
	size_t kpattern = bf.get_pattern(bf.pattern_hash(kmer));
	std::cerr << "kmer pattern number: " << kpattern << std::endl;
	std::cerr << "gibberish block: " << bf.get_block(bf.block_hash(gibberish)) << std::endl;
	std::cerr << "gibberish pattern number: " << bf.get_pattern(bf.pattern_hash(gibberish)) << std::endl;
	std::cerr << "kmerp: ";
	for(size_t i = kpattern; i < kpattern + 2;++i){
		for(size_t j = 0; j < 4; ++j){
			std::cerr << std::bitset<64>(bf.patterns.get()[i][j]) << ", ";	
		}
	}
	std::cerr << std::endl;

	std::cerr << "Table: \n";
	for(size_t i = 0; i < bf.size() / bits_per_char / 64; ++i){
		// if(i % (bf.block_size/bits_per_char) == 0){
		// 	std::cerr << "\nBLOCK" << std::endl;
		// }
		std::cerr << "block: ";
		for(size_t j = 0; j < 4; ++j){
			std::cerr << std::bitset<64>(bf.bit_table_.get()[2*i][j]) << ", ";
		}
		for(size_t j = 0; j < 4; ++j){
			std::cerr << std::bitset<64>(bf.bit_table_.get()[2*i+1][j]) << ", ";
		}
		std::cerr << std::endl;
	}
	std::cerr << std::endl;
	std::cerr << "bf contains kmer ? " << bf.contains(kmer) << std::endl;
	std::cerr << "bf contains gibberish ? " << bf.contains(gibberish) << std::endl;
	assert(bf.contains(kmer));


	return 0;
}