#include "readutils.hh"
#include <algorithm>
#include <iterator>
#include <iostream>

KSEQ_DECLARE(BGZF*)

namespace readutils{

	std::unordered_map<std::string, std::string> CReadData::rg_to_pu{};
	std::unordered_map<std::string, int> CReadData::rg_to_int{};

	CReadData::CReadData(bam1_t* bamrecord, bool use_oq){ //TODO: flag to use OQ
		this->name = bam_get_qname(bamrecord);
		this->seq = bam_seq_str(bamrecord);
		if(use_oq){
			const uint8_t* oqdata = bam_aux_get(bamrecord, "OQ"); // this will be null on error
			//we should throw in that case
			if(oqdata == NULL){
				std::cerr << "Error: --use-oq was specified but unable to read OQ tag " << 
				"on read " << this->name << std::endl;
				if(errno == ENOENT){
					std::cerr << "OQ not found. Try again without the --use-oq option." << std::endl;
				} else if(errno == EINVAL){
					std::cerr << "Tag data is corrupt. Repair the tags and try again." << std::endl;
				}
				throw std::invalid_argument("Unable to read OQ tag.");
			}
			const std::string oq = bam_aux2Z(oqdata);
			std::transform(oq.cbegin(), oq.cend(), std::back_inserter(this->qual),
				[](const char& c) -> uint8_t {return c - 33;});
		} else {
			std::copy(bam_get_qual(bamrecord), bam_get_qual(bamrecord) + bamrecord->core.l_qseq,
				std::back_inserter(this->qual));
		}
		if(bam_is_rev(bamrecord)){
			std::reverse(this->seq.begin(), this->seq.end());
			//can't do reversing and complementing at the same time because of how transform
			//works + the source/destination range overlap.
			std::transform(this->seq.begin(), this->seq.end(), this->seq.begin(),
				[](char c) -> char {return seq_nt16_str[seq_nt16_table[('0' + 3-seq_nt16_int[seq_nt16_table[c]])]];});
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
			throw std::invalid_argument("Unable to read RG tag.");;
		}
		this->rg = bam_aux2Z(bam_aux_get(bamrecord, "RG"));
		if(rg_to_pu.count(this->rg) == 0){
			rg_to_int[this->rg] = rg_to_int.size();
			// we just need something unique here.
			rg_to_pu[this->rg] = rg; //when loaded from the header this is actually a PU 
		}
		this->second = bamrecord->core.flag & BAM_FREAD2; // 0x80
		this->errors.resize(bamrecord->core.l_qseq);
	}
	// if second is >1, that means infer.

	CReadData::CReadData(kseq::kseq_t* fastqrecord, std::string rg, int second, std::string namedelimiter):
	seq(fastqrecord->seq.s), skips(seq.length(), false), errors(seq.length(), false)
	{
		// this->seq = std::string(fastqrecord->seq.s);
		// this->qual.assign(fastqrecord->qual.s, fastqrecord->qual.l);
		std::string quals(fastqrecord->qual.s);
		std::transform(quals.begin(), quals.end(), std::back_inserter(this->qual),
			[](char c) -> int {return c - 33;});
		// this->skips.resize(this->seq.length(), false);

		std::string fullname(fastqrecord->name.s);
		size_t current_pos = fullname.find(namedelimiter);
		std::string first_name = fullname.substr(0, current_pos);

		while(rg == "" && current_pos != std::string::npos){
			// if we need to find rg
			fullname = fullname.substr(current_pos+1); //get the right part, excluding the delimiter
			current_pos = fullname.find(namedelimiter); //reset the delimiter; this is npos if no more fields
			if(fullname.substr(0,3) == "RG:"){
				size_t last_colon = fullname.find_last_of(":", current_pos); // current_pos is last char to search
				rg = fullname.substr(last_colon+1, current_pos);
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
		// this->errors.resize(this->seq.length(), false);

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
		std::vector<bool> unskipped_errs;
		std::transform(this->skips.cbegin(), this->skips.cend(), this->errors.cbegin(),
			std::back_inserter(unskipped_errs),
			[](const bool &s, const bool &e) -> bool {return (!s) && e;});
		return unskipped_errs;
	}

	void CReadData::infer_read_errors(const bloom::Bloom& b, const std::vector<int>& thresholds, int k){
		std::array<std::vector<size_t>,2> overlapping = bloom::overlapping_kmers_in_bf(this->seq, b, k);
		std::vector<size_t> in = overlapping[0];
		std::vector<size_t> possible = overlapping[1];
		for(size_t i = 0; i < errors.size(); ++i){
			this->errors[i] = (in[i] <= thresholds[possible[i]] || this->qual[i] <= INFER_ERROR_BAD_QUAL);
#ifndef NDEBUG
			if(possible[i] > k){std::cerr << "seq:" << this->seq << " " << possible[i] << "WARNING: Invalid i: " << i << std::endl;}
			// if(i>=k-1 && std::string(this->seq, i-k+1, k) == "CCCCCCCCCTCGCCCCCCCCCCCCCCCCCCC"){
			// 	std::cerr << "seq: " << this->seq << "\n";
			// 	std::cerr << "in: ";
			// 	std::copy(in.begin()+i-k+1, in.begin()+i+1, std::ostream_iterator<size_t>(std::cerr, ", "));
			// 	std::cerr << "\npossible: ";
			// 	std::copy(possible.begin()+i-k+1, possible.begin()+i+1, std::ostream_iterator<size_t>(std::cerr, ", "));
			// 	std::cerr << "\nerrors: ";
			// 	std::copy(this->errors.begin()+i-k+1, this->errors.begin()+i+1, std::ostream_iterator<bool>(std::cerr, ", "));
			// 	std::cerr << std::endl;
			// }			
#endif
		}
	}

	size_t CReadData::correct_one(const bloom::Bloom& t, int k){
		int best_fix_len = 0;
		char best_fix_base;
		size_t best_fix_pos = std::string::npos;
		for(size_t i = 0; i < this->seq.length(); ++i){
			std::string original_seq(this->seq); //copy original
			for(const char& c : {'A','C','G','T'}){
				if(this->seq[i] == c){continue;}
				original_seq[i] = c;
				size_t start = i > k - 1 ? i - k + 1 : 0;
				//test a kmer to see whether its worth counting them all
				//i'm not sure any performance gain is worth it, but this is how Lighter does it
				size_t magic_start = i > k/2 - 1 ? std::min(i - k/2 + 1, original_seq.length()-k) : 0;
				bloom::Kmer magic_kmer(k);
				for(size_t j = magic_start; j <= magic_start + k - 1; ++j){
					magic_kmer.push_back(original_seq[j]);
				}
				//
				if(t.query(magic_kmer)){
					//94518
					int n_in = bloom::biggest_consecutive_trusted_block(
						original_seq.substr(start, 2*k - 1),t,k,best_fix_len);
#ifndef NDEBUG					
					std::cerr << "Found a kmer: " << magic_kmer << " i: " << i << " Fix len: " << n_in << std::endl;
#endif
					if(n_in > best_fix_len){ //94518
						best_fix_base = c;
						best_fix_pos = i;
						best_fix_len = n_in;
					} else if(n_in == best_fix_len && this->qual[i] < this->qual[best_fix_pos]){
						best_fix_base = c;
						best_fix_pos = i;
					}
				}
			}
		}
		if(best_fix_len > 0){
			this->seq[best_fix_pos] = best_fix_base;
		}
		return best_fix_pos;
	}

	//this is a chonky boi
	std::vector<bool> CReadData::get_errors(const bloom::Bloom& trusted, int k, int minqual, bool first_call){
		std::string original_seq(this->seq);
		size_t bad_prefix = 0;
		size_t bad_suffix = std::string::npos;
		bool multiple = false; //whether there were any ties
#ifndef NDEBUG
		std::cerr << "Correcting seq: " << original_seq << std::endl;
#endif
		std::array<size_t,2> anchor = bloom::find_longest_trusted_seq(this->seq, trusted, k);
#ifndef NDEBUG
		std::cerr << "Initial anchors: [" << anchor[0] << ", " << anchor[1] << "]" << std::endl;
#endif
		if(anchor[0] == std::string::npos){ //no trusted kmers in this read.
			multiple = true;
			size_t corrected_idx = this->correct_one(trusted, k);
			if(corrected_idx == std::string::npos){
				return this->errors;
			} else {
				anchor = bloom::find_longest_trusted_seq(this->seq, trusted, k);
#ifndef	NDEBUG
				std::cerr << "Created anchor: [" << anchor[0] << ", " << anchor[1] << "]" << std::endl;
#endif				
				this->errors[corrected_idx] = true;
			}
		}
		if(anchor[0] == 0 && anchor[1] == std::string::npos){ //all kmers are trusted
			return this->errors;
		}
		//we're guaranteed to have a valid anchor now.
		//min is in case anchor[1] is npos.
		size_t anchor_len = std::min(anchor[1], this->seq.length()-1) + 1 - anchor[0];
		bool corrected = false; //whether there were any corrections
		//right side
		if(anchor[1] != std::string::npos){
			//number of trusted kmers >= k
			if(anchor_len - k + 1 >= k){ //number of trusted kmers 
				bool current_multiple;
				std::tie(anchor[1], current_multiple) = bloom::adjust_right_anchor(anchor[1], this->seq, trusted, k);
#ifndef NDEBUG
				std::cerr << "Adjust R Multiple: " << current_multiple << std::endl;
#endif
				multiple = multiple || current_multiple;
			}
			for(size_t i = anchor[1] + 1; i < this->seq.length();){
				size_t start = i - k + 1; //seq containing all kmers that are affected
				std::vector<char> fix;
				size_t fixlen;
				bool current_multiple;
				std::tie(fix, fixlen, current_multiple) = bloom::find_longest_fix(this->seq.substr(start, std::string::npos), trusted, k);
#ifndef NDEBUG
				std::cerr << "R fix Multiple: " << current_multiple << std::endl;
#endif
				multiple = multiple || current_multiple;
				size_t next_untrusted_idx = start + fixlen;
				if(next_untrusted_idx > i){
					if(fix.size() > 1){
						// multiple = true; //we take care of this in function
						// i = index of erroneous base
						//to = ( i + kmerLength - 1 < readLength ) ? i + kmerLength - 1 : readLength - 1 ;
						//to = inclusive biggest index the fix can possibly go to. (= largest_possible_idx)
						//maxto = largest index the fix actually goes to INCLUSIVE; next_untrusted_idx-1
						//if( maxTo <= to || to - i + 1 < kmerLength ) ...
						//i + k
						multiple = true;
						size_t largest_possible_idx = std::min(i + k - 1, this->seq.length()-1);
#ifndef NDEBUG
						std::cerr << "next_untrusted_idx: " << next_untrusted_idx << std::endl;
						std::cerr << "largest_possible_idx (to): " << largest_possible_idx << std::endl;
						std::cerr << "largest_possible_idx-i+1: " << largest_possible_idx-i+1 << std::endl;
#endif					//20629746
						if(next_untrusted_idx-1 <= largest_possible_idx || largest_possible_idx - i + 1 < k){
							// size_t trimstart = next_untrusted_idx;
							bad_suffix = i; //readlength - trimstart
#ifndef NDEBUG
							//if there's a tie and we haven't gone the max number of kmers, end correction
							std::cerr << "i: " << i << " next_untrusted_idx " << next_untrusted_idx << std::endl;
							std::cerr << "Fixlen is " << fixlen << " Fix: ";
							for(char c : fix){
								std::cerr << c << " ";
							}
							std::cerr << std::endl;
							std::cerr << "Tie and fix not long enough; ending correction early!" << std::endl;
#endif
							break;
						}
					} else {
						this->seq[i] = fix[0];
						this->errors[i] = true;
					}
					corrected = true;
#ifndef NDEBUG
					std::cerr << "Error detected at position " << i << ". Advancing " << fixlen - k + 1 << "." << std::endl;
#endif
					i += fixlen - k + 1; // i = next_untrusted_idx
				} else {
					//couldn't find a fix; skip ahead and try again if long enough
#ifndef NDEBUG
					std::cerr << "Couldn't fix position " << i << " Cutting read and trying again." << std::endl;//". Skipping ahead " << k - 1 << "." << std::endl;
#endif
					bad_suffix = i;
					break;
					// i += k-1; //move ahead and make i = k-1 as the first base in the new kmer
					// if(this->seq.length() - i + k <= (seq.length()/2) || this->seq.length() - i + k <= 2*k ){
					// 	//sequence not long enough. end this side.
					// 	break;
					// }
				}
			}
		}
		//left side
		if(anchor[0] != 0){
			//the bad base is at anchor[0]-1, then include the full kmer for that base.
			// std::string sub = this->seq.substr(0, anchor[0] - 1 + k);
			std::string revcomped(this->seq.length(), 'N');
			std::transform(this->seq.rbegin(), this->seq.rend(), revcomped.begin(),
				[](char c) -> char {return seq_nt16_str[seq_nt16_table[('0' + 3-seq_nt16_int[seq_nt16_table[c]])]];});
			//if num of trusted kmers >= k, see if anchor needs adjusting.
			if(anchor_len - k + 1 >= k){ 
				size_t left_adjust;
				bool current_multiple;
				std::tie(left_adjust, current_multiple) = bloom::adjust_right_anchor(revcomped.length()-anchor[0]-1, revcomped, trusted, k);
				anchor[0] = revcomped.length()-left_adjust-1; //change back to original coordinates
#ifndef NDEBUG
				std::cerr << "Adjust L Multiple: " << current_multiple << std::endl;
#endif
				multiple = multiple || current_multiple;
			}
			//
			for(int i = anchor[0] - 1; i >= 0;){ //index of erroneous base in original seq
				int j = revcomped.length()-i-1; //index of erroneous base in reversed seq
				size_t start = j - k + 1; //seq containing all kmers that are affected
				//but [j -k + 1, npos) in reverse space.
				std::string sub = revcomped.substr(start, std::string::npos); //get the right subsequence
				std::vector<char> fix;
				size_t fixlen;
				bool current_multiple;
				std::tie(fix, fixlen, current_multiple) = bloom::find_longest_fix(sub, trusted, k, true); //155392 TODO: add reverse_test to adjust_right_anchor
#ifndef NDEBUG
				std::cerr << "L Fix Multiple: " << current_multiple << std::endl;
#endif
				multiple = multiple || current_multiple;
				//new i is return value + start // i += r -k + 1
				//new j should be return value + start // j += r - k + 1
				// new j should be fixlen + start; j += fixlen - k + 1
				// j = fixlen + start; j = j - k + 1 + fixlen; j += fixlen -k + 1
				size_t next_untrusted_idx = start + fixlen; // next untrusted idx in j space
				if( next_untrusted_idx > j){
					if(fix.size() > 1){
						multiple = true; // we don't have to do this here because we do it in function
						//this value needs to be fixed i think
						//to = ( i - kmerLength + 1 < 0 ) ? 0 : ( i - kmerLength + 1 ) ; 
						//( minTo >= to || i - to + 1 < kmerLength )
						//Line num: 21537
						// if(next_untrusted_idx - 1 <= std::max(start + (size_t)2*k - 1, revcomped.length())){
						size_t largest_possible_idx = std::min(j + (size_t)k - 1, revcomped.length()-1);
						if(next_untrusted_idx-1 <= largest_possible_idx || largest_possible_idx - j + 1 < k){ //595573
							//if there's a tie and we haven't gone the max number of kmers, end correction
							bad_prefix = i;
							break;
						}
					} else {
						revcomped[j] = fix[0];
						this->errors[i] = true;
					}
					corrected = true;
#ifndef NDEBUG
					std::cerr << "next_untrusted_idx: " << next_untrusted_idx << " j: " << j << std::endl;
					std::cerr << "Error detected at position " << i << ". Advancing " << next_untrusted_idx-j << " " << (fixlen - k + 1) << "." << std::endl;
#endif					
					i -= next_untrusted_idx - j; //fixlen - k + 1;
				} else {
					//couldn't find a fix; skip ahead and try again if long enough
#ifndef NDEBUG
					std::cerr << "Couldn't fix position " << i << "Cutting read and trying again." << std::endl;//". Skipping ahead " << k - 1 << "." << std::endl;
#endif			
					bad_prefix = i; //the last position in the new read
					break;
					// i -= k-1;
					// if(i + k <= (seq.length()/2) || i + k <= 2*k ){
					// 	//sequence not long enough. end this side.
					// 	break;
					// }
				}
			}
		}
#ifndef NDEBUG
		std::cerr << "Anchors: [" << anchor[0] << ", " << anchor[1] << "]" << std::endl; 
#endif


		// check for overcorrection and fix it
		if(corrected){
			bool adjust = true;
			//check that no trusted kmers were "fixed"
			bloom::Kmer kmer(k);
			size_t trusted_start = std::string::npos;
			size_t trusted_end = std::string::npos;
			for(size_t i = 0; i < original_seq.length() && adjust == true; ++i){
				kmer.push_back(original_seq[i]);
				if(kmer.valid() && trusted.query(kmer)){
					trusted_start = std::min(trusted_start,i-k+1);
					trusted_end = i;
					// std::cerr << "Trusted: " << trusted_start << " " << trusted_end << std::endl;
				} else {
					if(i > trusted_end){
						//we clear everything from beginning to end
						for(size_t j = trusted_start; j <= trusted_end; ++j){
							if(this->errors[j]){
								adjust = false;
								break;
							}
						}
						trusted_start = std::string::npos;
						trusted_end = std::string::npos;
					}
				}
			}
			//if we get to the end but have a trusted block
			//we don't actually trust that block :|
			// if(trusted_end != std::string::npos){
			// 	for(size_t j = trusted_start; j <= trusted_end; ++j){
			// 		if(this->errors[j]){
			// 			std::cerr << "(2) Problem: j: " << j << " " << trusted_start << " " << trusted_end << std::endl;
			// 			adjust = false;
			// 			break;
			// 		}
			// 	}
			// }
#ifndef NDEBUG
			std::cerr << "Adjust: " << adjust << " Multiple: " << multiple << std::endl;
#endif
			adjust = adjust && !multiple;//made it through the loop and no ties during correction
			// std::cerr << "Read corrected. Adjust threshold? " << adjust << std::endl;
#ifndef NDEBUG
			std::cerr << "Errors before adjustment: ";
			for(const bool& b: this->errors){
				std::cerr << b;
			}
			std::cerr << std::endl;
			std::cerr << "Seq After Correction: " << seq << std::endl;
#endif
			int ocwindow = 20;
			int base_threshold = 4;
			int threshold = base_threshold;
			double occount = 0;
			//check for overcorrection
			std::vector<int> overcorrected_idx; //push_back overcorrected indices in order
			//then from overcorrected.begin() - k to overcorrected.end() + k should all be reset.
			for(int i = 0; i < this->seq.length(); ++i){
				if(this->errors[i] && seq_nt16_int[seq_nt16_table[original_seq[i]]] < 4){ //increment correction count
					if(this->qual[i] <= minqual){
						occount += 0.5;
					} else {
						++occount;
					}
				}
				if(i >= ocwindow && this->errors[i-ocwindow] && seq_nt16_int[seq_nt16_table[original_seq[i-ocwindow]]] < 4){ //decrement count for not in window
					if(this->qual[i-ocwindow] <= minqual){
						occount -= 0.5;
					} else {
						--occount;
					}
				}
				//set threshold
				threshold = adjust && i >= ocwindow && i + ocwindow - 1 < this->seq.length() ? 
					base_threshold + 1 : base_threshold;
				//determine if overcorrected
#ifndef NDEBUG
				std::cerr << "Occount: " << occount << " Threshold: " << threshold;
				std::cerr << " Seq: " << original_seq[i] << " (" << seq_nt16_int[seq_nt16_table[original_seq[i]]];
				std::cerr << ")" << " Q: " << +this->qual[i] << std::endl;
#endif
				if(occount > threshold && this->errors[i]){
					overcorrected_idx.push_back(i);
				}

			}
#ifndef NDEBUG
			std::cerr << "Overcorrected indices (" << overcorrected_idx.size() << "): ";
			std::copy(overcorrected_idx.begin(), overcorrected_idx.end(), std::ostream_iterator<int>(std::cerr, ", "));
			std::cerr << std::endl;
#endif
			//Line num: 23026
			for(int oc_idx : overcorrected_idx){
				if(this->errors[oc_idx]){ //overcorrected idx hasn't been addressed yet
					int start = oc_idx-k+1; //the beginningmost position to check
					start = start >= 0? start : 0;
					int end = oc_idx+k; //the endmost position to check
					end = end < this->seq.length() ? end : this->seq.length();
					//we start iteration but we need to unfix anything within k of an overcorrected position
					//OR within k of one of those fixed positions.
					for(int i = start; i < end; ++i){
						if(this->errors[i]){
							this->errors[i] = false;
							//changint the end must come before the start change because the start
							//change changes i!
							if(i+k > end){ //change the end if we need to
								end = i+k < this->seq.length() ? i+k : this->seq.length();
							}
							//i will be 1 greater than start, so rather than i-k+1 we have i-k.
							if(i-k < start){ //go back a bit if we need to; +1 comes from the loop
								i = i-k+1 >= 0 ? i-k : -1; //+1 will come from the loop
								start = i;
							}
						}
					}
				}
			}
		}
		if(first_call && bad_prefix > 0 && (bad_prefix >= this->seq.length() / 2 || bad_prefix >= 2*k)){
#ifndef NDEBUG			
			std::cerr << "bad_prefix: " << bad_prefix << std::endl;
#endif
			CReadData subread = this->substr(0, bad_prefix+1); //2nd argument is length
			std::vector<bool> suberrors = subread.get_errors(trusted, k, minqual, false);
			std::copy(suberrors.begin(), suberrors.end(), this->errors.begin());
		}
		if(first_call && bad_suffix < std::string::npos && bad_suffix < this->seq.length() &&
		(this->seq.length()-bad_suffix > this->seq.length()/2 || this->seq.length()-bad_suffix > 2*k)){
#ifndef NDEBUG
			std::cerr << "bad_suffix: " << bad_suffix << std::endl;
#endif
			CReadData subread = this->substr(bad_suffix, std::string::npos);
			std::vector<bool> suberrors = subread.get_errors(trusted, k, minqual, false);
			std::copy(suberrors.begin(), suberrors.end(), this->errors.begin()+bad_suffix);
		}
		if(std::find(this->errors.begin(), this->errors.end(), true) != this->errors.end()){
			corrected = true;
		}

		this->seq = original_seq;
		return this->errors;
	}

	std::vector<uint8_t> CReadData::recalibrate(const covariateutils::dq_t& dqs, int minqual) const{
		//TODO this function uses a narrowing conversion which is unsafe
		std::vector<uint8_t> recalibrated(this->qual);
		int rg = this->get_rg_int();
		for(int i = 0; i < this->seq.length(); ++i){
			uint8_t q = this->qual[i];
			if(q >= minqual){
				recalibrated[i] = dqs.meanq[rg] + dqs.rgdq[rg] + dqs.qscoredq[rg][q] +
					dqs.cycledq[rg][q][this->second][i];
				if(i > 0){
					int first = seq_nt16_int[seq_nt16_table[this->seq[i-1]]];
					int second = seq_nt16_int[seq_nt16_table[this->seq[i]]];
					if(first < 4 && second < 4){
						int8_t dinuc = covariateutils::dinuc_to_int(first, second);
						recalibrated[i] += dqs.dinucdq[rg][q][dinuc];
					}
				}
			}
		}
		return recalibrated;
	}

	CReadData CReadData::substr(size_t pos, size_t count) const{
		CReadData ret = (*this);
		ret.seq = ret.seq.substr(pos, count);
		size_t len = ret.seq.size();
		ret.qual.resize(len);
		ret.skips.resize(len);
		ret.errors.resize(len);
		std::copy(this->qual.cbegin()+pos, this->qual.cbegin()+pos+len, ret.qual.begin());
		std::copy(this->skips.cbegin()+pos, this->skips.cbegin()+pos+len, ret.skips.begin());
		std::copy(this->errors.cbegin()+pos, this->errors.cbegin()+pos+len, ret.errors.begin());
		return ret;
	}


}

