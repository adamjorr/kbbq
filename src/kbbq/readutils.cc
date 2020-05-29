#include "readutils.hh"
#include <algorithm>
#include <iostream>

KSEQ_DECLARE(BGZF*)

namespace readutils{

	std::unordered_map<std::string, std::string> CReadData::rg_to_pu{};
	std::unordered_map<std::string, int> CReadData::rg_to_int{};

	CReadData::CReadData(bam1_t* bamrecord, bool use_oq){ //TODO: flag to use OQ
		this->name = bam_get_qname(bamrecord);
		this->seq = bam_seq_str(bamrecord);
		if(use_oq){
			uint8_t* oqdata = bam_aux_get(bamrecord, "OQ"); // this will be null on error
			//we should throw in that case
			if(oqdata == NULL){
				std::cerr << "Error: --use-oq was specified but unable to read OQ tag " << 
				"on read " << this->name << std::endl;
				if(errno == ENOENT){
					std::cerr << "OQ not found. Try again without the --use-oq option." << std::endl;
				} else if(errno == EINVAL){
					std::cerr << "Tag data is corrupt. Repair the tags and try again." << std::endl;
				}
				throw;
			}
			std::string oq(bam_aux2Z(oqdata));
			std::transform(oq.begin(), oq.end(), std::back_inserter(this->qual),
				[](char c) -> int {return c - 33;});
		} else {
			this->qual.assign(bam_get_qual(bamrecord), bam_get_qual(bamrecord) + bamrecord->core.l_qseq);
		}
		if(bam_is_rev(bamrecord)){
			std::reverse(this->seq.begin(), this->seq.end());
			std::reverse(this->qual.begin(), this->qual.end());
		}
		this->skips.resize(bamrecord->core.l_qseq, 0);
		uint8_t* rgdata = bam_aux_get(bamrecord, "RG");
		if(rgdata == NULL){
			std::cerr << "Error: Unable to read RG tag on read " << this->name << std::endl;
			if(errno == ENOENT){
				std::cerr << "RG not found. " <<
				"Every read in the BAM must have an RG tag; add tags with " <<
				"samtools addreplacerg and try again." << std::endl;
			} else if(errno == EINVAL){
				std::cerr << "Tag data is corrupt. Repair the tags and try again." << std::endl;
			}
			throw;
		}
		this->rg = bam_aux2Z(bam_aux_get(bamrecord, "RG"));
		if(rg_to_pu.count(this->rg) == 0){
			rg_to_int[this->rg] = rg_to_int.size();
			// we just need something unique here.
			rg_to_pu[this->rg] = rg; //when loaded from the header this is actually a PU 
		}
		this->second = bamrecord->core.flag & 0x80;
		this->errors.resize(bamrecord->core.l_qseq);
	}
	// if second is >1, that means infer.

	CReadData::CReadData(kseq::kseq_t* fastqrecord, std::string rg, int second, std::string namedelimiter){
		this->seq = std::string(fastqrecord->seq.s);
		// this->qual.assign(fastqrecord->qual.s, fastqrecord->qual.l);
		std::string quals(fastqrecord->qual.s);
		std::transform(quals.begin(), quals.end(), std::back_inserter(this->qual),
			[](char c) -> int {return c - 33;});
		this->skips.resize(this->seq.length(), false);

		std::string fullname(fastqrecord->name.s);
		size_t current_pos = fullname.find(namedelimiter);
		std::string first_name = fullname.substr(0, current_pos);

		while(rg == "" && current_pos != std::string::npos){
			// if we need to find rg
			fullname = fullname.substr(current_pos+1); //get the right part, excluding the delimiter
			current_pos = fullname.find(namedelimiter); //reset the delimiter; this is npos if no more fields
			if(fullname.substr(0,3) == "RG:"){
				size_t last_colon = fullname.find_last_of(":", current_pos); // current_pos is last char to search
				rg = fullname.substr(last_colon, current_pos);
			}
		}
		this->rg = rg;

		if(second > 1){
			std::string tail = first_name.substr(first_name.length() - 2);
			second = (tail == "/2");
			if(second || tail == "/1"){
				first_name = first_name.substr(0, first_name.length() - 2);
			}
		}
		this->name = first_name;
		this->second = second;
		this->errors.resize(this->seq.length(), false);

		if(rg_to_pu.count(this->rg) == 0){
			rg_to_int[this->rg] = rg_to_int.size();
			rg_to_pu[this->rg] = rg; //when loaded from the header this is actually a PU.
		}
	}

	void CReadData::load_rgs_from_bamfile(bam_hdr_t* header){
		std::string hdrtxt(header->text);
		size_t linedelim = hdrtxt.find('\n');
		// kstring_t val;
		// while (sam_hdr_find_tag_pos(header, "RG", i, "ID", val) == 0){
		// 	std::string rgid(val.s);
		// 	sam_hdr_find_tag_pos(header, "RG", i, "PU", val);
		// 	std::string pu(val.s);
		// 	if(rg_to_pu.count(rgid) == 0){
		// 		rg_to_int[rgid] = rg_to_int.size();
		// 		rg_to_pu[rgid] = pu;
		// 	}
		// }
		while(linedelim != std::string::npos){
			if(hdrtxt.substr(0,3) == "@RG"){
				std::string id("");
				std::string pu("");
				std::string line = hdrtxt.substr(0, linedelim);
				size_t tokendelim = line.find_first_of("\t ");
				while(tokendelim != std::string::npos){
					if(line.substr(0,3) == "ID:"){
						//id
						id = line.substr(4, tokendelim);
					} else if (line.substr(0, 3) == "PU:"){
						//pu
						pu = line.substr(4, tokendelim);
					}
					line = line.substr(tokendelim);
					tokendelim = line.find_first_of("\t ");
				}
				if(id != ""){
					rg_to_int[id] = rg_to_int.size();
					rg_to_pu[id] = pu;
				}
			}
			hdrtxt = hdrtxt.substr(linedelim);
			linedelim = hdrtxt.find('\n');
		}
	}

	std::string CReadData::str_qual(){
		std::string str_qual;
		for(size_t i = 0; i < this->qual.size(); i++){
			str_qual.push_back(this->qual[i] + 33);
		}
		return str_qual;
	}

	std::string CReadData::canonical_name(){
		std::string suffix;
		if (this->second){
			suffix = "/2";
		}
		else{
			suffix = "/1";
		}
		return this->name + suffix;
	}

	std::vector<bool> CReadData::not_skipped_errors() const{
		std::vector<bool> unskipped = this->skips;
		unskipped.flip();
		for(size_t i = 0; i < unskipped.size(); i++){
			unskipped[i] = (unskipped[i] && this->errors[i]);
		}
		return unskipped;
	}

	void CReadData::infer_read_errors(const bloom::bloomary_t& b, const std::vector<int>& thresholds, int k){
		std::array<std::vector<size_t>,2> overlapping = bloom::overlapping_kmers_in_bf(this->seq, b, k);
		std::vector<size_t> in = overlapping[0];
		std::vector<size_t> possible = overlapping[1];
		for(size_t i = 0; i < errors.size(); ++i){
			// std::cerr << "in " << i << ": " << in[i] << std::endl;
			// std::cerr << "possible " << i << ": " << possible[i] << std::endl;
			if(possible[i] > k){std::cerr << "seq:" << this->seq << " " << possible[i] << "WARNING: Invalid i: " << i << std::endl;}
			this->errors[i] = (in[i] <= thresholds[possible[i]] || this->qual[i] <= 6);
		}
	}

	int CReadData::correct_one(const bloom::bloomary_t& t, int k){
		int best_fix_len = 0;
		char best_fix_base;
		size_t best_fix_pos = 0;
		for(size_t i = 0; i < this->seq.length(); ++i){
			std::string original_seq(this->seq); //copy original
			for(const char* c : {"A","C","G","T"}){
				original_seq[i] = *c;
				size_t start_pos = (size_t) std::min(std::max((int)i-k/2+1,0), (int)this->seq.length()-k);
				int n_in = bloom::nkmers_in_bf(this->seq.substr(start_pos),t,k);
				if(n_in > best_fix_len){
					best_fix_base = *c;
					best_fix_pos = i;
					best_fix_len = n_in;
				}
			}
		}
		if(best_fix_len > 0){
			this->seq[best_fix_pos] = best_fix_base;
		}
		return best_fix_len;
	}

	//this is a chonky boi
	void CReadData::get_errors(const bloom::bloomary_t& trusted, int k, int minqual){
		std::string original_seq(this->seq);
		std::array<size_t,2> anchor = bloom::find_longest_trusted_seq(this->seq, trusted, k);
		if(anchor[0] == std::string::npos){ //no trusted kmers in this read.
			if(this->correct_one(trusted, k) == 0){
				return;
			} else {
				anchor = bloom::find_longest_trusted_seq(this->seq, trusted, k);
			}
		}
		if(anchor[0] == 0 && anchor[1] == std::string::npos){ //all kmers are trusted
			return;
		}
		//we're guaranteed to have a valid anchor now.
		bool corrected = false; //whether there were any corrections
		bool multiple = false; //whether there were any ties
		//right side
		if(anchor[1] != std::string::npos){
			for(size_t i = anchor[1] + 1; i < this->seq.length();){
				size_t start = i - k + 1; //seq containing all kmers that are affected
				std::pair<std::vector<char>, int> lf = bloom::find_longest_fix(this->seq.substr(start, std::string::npos), trusted, k);
				std::vector<char> fix = std::get<0>(lf);
				int fixlen = std::get<1>(lf); //new i is return value + start // i += r -k + 1
				if(fixlen-k+1 > 0){
					this->seq[i] = fix[0];
					this->errors[i] = true;
					i += fixlen - k + 1;
					corrected = true;
					if(fix.size() > 1){
						multiple = true;
					}
				} else {
					//couldn't find a fix; skip ahead and try again if long enough
					i += k-1; //move ahead and make i = k-1 as the first base in the new kmer
					if(this->seq.length() - i + k <= (seq.length()/2) || this->seq.length() - i + k <= 2*k ){
						//sequence not long enough. end this side.
						break;
					}
				}
			}
		}
		//left side
		if(anchor[0] != 0){
			//the bad base is at anchor[0]-1, then include the full kmer for that base.
			std::string sub = this->seq.substr(0, anchor[0] - 1 + k);
			std::string revcomped;
			for(auto it = sub.rbegin(); it != sub.rend(); it++){
					int c = seq_nt4_table[*it];
					char d = c < 4 ? seq_nt16_str[seq_nt16_table[3 - c]] : c;
					revcomped.push_back(d); //complement then turn to 4-bit encoded then turn to char
				}
			for(int i = anchor[0] - 1; i > 0;){ //iterating in original seq space
				int j = anchor[0] - 1 + k - 1 - i; //iterating in reversed space
				int end = i + k; //seq containing all kmers that are affected is [0, end) in original space
				//but [j -k + 1, npos) in reverse space. 
				std::string sub = revcomped.substr(j - k + 1, std::string::npos); //get the right subsequence
				std::pair<std::vector<char>, int> lf = bloom::find_longest_fix(sub, trusted, k);
				std::vector<char> fix = std::get<0>(lf);
				int fixlen = std::get<1>(lf); //new i is return value + start // i += r -k + 1
				//new j should be return value + start // j += r - k + 1
				if(fixlen-k+1 > 0){
					int c = seq_nt4_table[fix[0]];
					char d = c < 4 ? seq_nt16_str[seq_nt16_table[3 - c]] : c;
					revcomped[j] = d;
					this->errors[i] = true;
					i -= fixlen - k + 1;
					corrected = true;
					if(fix.size() > 1){
						multiple = true;
					}
				} else {
					//couldn't find a fix; skip ahead and try again if long enough
					i -= k+1;
					if(i + k <= (seq.length()/2) || i + k <= 2*k ){
						//sequence not long enough. end this side.
						break;
					}
				}
			}
		}
		// check for overcorrection and fix it
		if(corrected){
			bool adjust = !multiple; //if there were any ties, don't adjust
			int ocwindow = 20;
			int threshold = 4;
			double occount = 0;
			//check for overcorrection
			std::vector<int> overcorrected_idx; //push_back overcorrected indices in order
			//then from overcorrected.begin() - k to overcorrected.end() + k should all be reset.
			for(int i = 0; i < this->seq.length(); ++i){
				if(this->errors[i] && seq_nt4_table[this->seq[i]] < 4){ //increment correction count
					if(this->qual[i] <= minqual){
						occount += 0.5;
					} else {
						++occount;
					}
				}
				if(i > ocwindow && this->errors[i-ocwindow] && seq_nt4_table[this->seq[i-ocwindow]] < 4){ //decrement count for not in window
					if(this->qual[i-ocwindow] <= minqual){
						occount -= 0.5;
					} else {
						--occount;
					}
				}
				//set threshold
				if(adjust && i >= ocwindow){
					threshold++;
				}
				if(adjust && i + ocwindow - 1 < this->seq.length()){
					threshold--;
				}
				//determine if overcorrected
				if(occount > threshold && this->errors[i]){
					overcorrected_idx.push_back(i);
				}

			}

			if(overcorrected_idx.size() > 0){
				int start = overcorrected_idx[0]-k+1; //the beginningmost position to check
				start = start >= 0? start : 0;
				int end = overcorrected_idx.back()+k-1; //the endmost position to check
				end = end < this->seq.length() ? end : this->seq.length();
				//we start iteration but we need to unfix anything within k of an overcorrected position
				//OR within k of one of those fixed positions.
				for(int i = start; i < end; ++i){
					if(this->errors[i]){
						this->errors[i] = false;
						if(i-k+1 < start){ //go back a bit if we need to
							i = i-k+1 >= 0 ? i-k+1 : 0;
							start = i;
						}
						if(i+k-1 > end){ //change the end if we need to
							end = i+k-1 < this->seq.length() ? i+k-1 : this->seq.length();
						}
					}
				}
			}
		}
		this->seq = original_seq;
	}

	std::vector<int> CReadData::recalibrate(const covariateutils::dq_t& dqs, int minqual) const{
		std::vector<int> recalibrated(this->qual);
		int rg = this->get_rg_int();
		for(int i = 0; i < this->seq.length(); ++i){
			int q = this->qual[i];
			if(q >= minqual){
				recalibrated[i] = dqs.meanq[rg] + dqs.rgdq[rg] + dqs.qscoredq[rg][q] +
				dqs.cycledq[rg][q][this->second][i];
				if(i > 0){
					int first = seq_nt4_table[this->seq[i-1]];
					int second = seq_nt4_table[this->seq[i]];
					if(first < 4 && second < 4){
						int8_t dinuc = covariateutils::dinuc_to_int(first, second);
						recalibrated[i] += dqs.dinucdq[rg][q][dinuc];
					}
				}
			}
		}
		return recalibrated;
	}
}

