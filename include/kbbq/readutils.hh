#ifndef READUTILS_H
#define READUTILS_H

#include <vector>
#include <unordered_map>
#include <string>
#include <errno.h>
//
#include <htslib/sam.h>
#include <htslib/kseq.h>
#include <htslib/bgzf.h>
#include "bloom.hh"
#include "covariateutils.hh"
#include "kseq.hh"

//fwd declare
namespace covariateutils{
	struct dq_t;
}

namespace readutils{

	static int8_t complement[16] = { 0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15 };

	// int correction_len(const bloom::bloomary_t& t, int k);

	//get fwd-orientation read string
	inline std::string bam_seq_str(bam1_t* bamrecord){
		std::string seq;
		unsigned char* s = bam_get_seq(bamrecord);
		for(int i = 0; i < bamrecord->core.l_qseq; ++i){
			seq.push_back(bam_is_rev(bamrecord) ? seq_nt16_str[complement[bam_seqi(s, i)]] : seq_nt16_str[bam_seqi(s, i)]);
		}
		return seq;
	}

	class CReadData{
		private:
			static std::unordered_map<std::string, std::string> rg_to_pu;
			static std::unordered_map<std::string, int> rg_to_int;
		public:
			CReadData(){}
			CReadData(bam1_t* bamrecord, bool use_oq = false);
			// hello?
			CReadData(kseq::kseq_t* fastqrecord, std::string rg = "", int second = 2, std::string namedelimiter = "_");
			std::string seq;
			std::vector<int> qual;
			std::vector<bool> skips;
			std::string name;
			std::string rg;
			bool second;
			std::vector<bool> errors;
			std::string str_qual();
			std::string canonical_name();
			inline int get_rg_int() const{return this->rg_to_int[this->rg];}
			inline std::string get_pu() const{return this->rg_to_pu[this->rg];}
			std::vector<bool> not_skipped_errors() const;
			//fill errors attribute given the bloom filter and thresholds.
			void infer_read_errors(const bloom::bloomary_t& b, const std::vector<int>& thresholds, int k);
			//fix one error and return the length of kmers it fixed
			int correct_one(const bloom::bloomary_t& t, int k);
			static void load_rgs_from_bamfile(bam_hdr_t* header);
			void get_errors(const bloom::bloomary_t& trusted, int k, int minqual = 6);
			std::vector<int> recalibrate(const covariateutils::dq_t& dqs, int minqual = 6) const;

	};
}

#endif
