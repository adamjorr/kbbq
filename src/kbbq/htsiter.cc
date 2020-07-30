#include "htsiter.hh"

namespace htsiter{

int BamFile::next(){return sam_itr_next(sf, itr, r);}
	// return next read as a string. if there are no more, return the empty string.
std::string BamFile::next_str(){return this->next() >= 0 ? readutils::bam_seq_str(r) : "";}
//
readutils::CReadData BamFile::get(){return readutils::CReadData(this->r, use_oq);}
//
void BamFile::recalibrate(const std::vector<uint8_t>& qual){
	uint8_t* q = bam_get_qual(this->r);
	if(set_oq){
		std::string qstr;
		std::transform(q, q + this->r->core.l_qseq, std::back_inserter(qstr),
			[](uint8_t c) -> char {return c + 33;}); //qual value to actual str
		//returns 0 on success, -1 on fail. We should consider throwing if it fails.
		if(bam_aux_update_str(this->r, "OQ", qstr.length()+1, qstr.c_str()) != 0){
			if(errno == ENOMEM){
				std::cerr << "Insufficient memory to expand bam record." << std::endl;
			} else if(errno == EINVAL){
				std::cerr << "Tag data is corrupt. Repair the tags and try again." << std::endl;
			}
			throw std::invalid_argument("Unable to update OQ tag.");
		}
	}
	if(bam_is_rev(this->r)){
		std::reverse_copy(qual.begin(), qual.end(), q);
	} else {
		std::copy(qual.begin(), qual.end(), q);
	}
	// for(int i = 0; i < this->r->core.l_qseq; ++i){
	// 	q[i] = (char)qual[i];
	// }
}
// TODO:: add a PG tag to the header
int BamFile::open_out(std::string filename){this->of = sam_open(filename.c_str(), "wb"); return sam_hdr_write(this->of, this->h);}
//
int BamFile::write(){return sam_write1(this->of, this->h, this->r);}

// FastqFile class

int FastqFile::next(){
	return kseq_read(r);
}

std::string FastqFile::next_str(){
	return this->next() >= 0? std::string(this->r->seq.s): "";
}

readutils::CReadData FastqFile::get(){
	return readutils::CReadData(this->r);
}

void FastqFile::recalibrate(const std::vector<uint8_t>& qual){
	for(int i = 0; i < this->r->qual.l; ++i){
		this->r->qual.s[i] = (char)(qual[i]+33);
	}
}

int FastqFile::open_out(std::string filename){
	ofh = bgzf_open(filename.c_str(),"w"); //mode should be "wu" if uncompressed output is desired.
	return ofh == 0 ? -1 : 0;
}

int FastqFile::write(){
	std::string name = ks_c_str(&this->r->name);
	std::string seq = ks_c_str(&this->r->seq);
	std::string comment = ks_c_str(&this->r->comment);
	std::string qual = ks_c_str(&this->r->qual);
	// std::string name = this->r->name.s ? std::string(this->r->name.s) : "";
	// std::string seq = this->r->seq.s ? std::string(this->r->seq.s) : "";
	// std::string comment = this->r->comment.s ? std::string(this->r->comment.s) : "";
	// std::string qual = this->r->qual.s ? std::string(this->r->qual.s) : "";
	std::string s("@" + name + "\n" + seq + "\n+" + comment + "\n" + qual + "\n");
	return bgzf_write(ofh, s.c_str(), s.length());
}

// KmerSubsampler
bloom::Kmer KmerSubsampler::next_kmer(){
	if(cur_kmer < kmers.size()){
		return kmers[cur_kmer++]; //return the current kmer and advance
	} else {
		readseq = file->next_str();
		kmer.reset();
		if(readseq.empty()){
			this->not_eof = false;
			return kmer; //no more sequences
		} else {
			kmers.clear();
			for(size_t i = 0; i < readseq.length(); ++i){
				kmer.push_back(readseq[i]);
				if(i >= k-1){
					kmers.push_back(kmer);
				}
			}
			cur_kmer = 0; //reset current kmer
			total_kmers += kmers.size();
			return this->next_kmer(); //try again
		}
	}
}

bloom::Kmer KmerSubsampler::next(){
	bloom::Kmer kmer = this->next_kmer();
	if(this->not_eof){
#ifdef KBBQ_USE_RAND_SAMPLER
		double p = std::rand() / (double)RAND_MAX;
		if(p < this->d.p()){
#else
		if(d(rng)){ //sampled
#endif
			return kmer;
		}
		else{ //try again
			return this->next();
		}
	}
	return kmer; // empty
}








}
