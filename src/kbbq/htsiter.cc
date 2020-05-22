#include "htsiter.hh"

namespace htsiter{

int BamFile::next(){return sam_itr_next(sf, itr, r);}
	// return next read as a string. if there are no more, return the empty string.
std::string BamFile::next_str(){return this->next() >= 0 ? readutils::bam_seq_str(r) : "";}
//
readutils::CReadData BamFile::get(){return readutils::CReadData(this->r, use_oq);}
//
void BamFile::recalibrate(std::vector<int> qual){
	char* q = (char *)bam_get_qual(this->r);
	if(set_oq){
		std::string qstr(q, this->r->core.l_qseq);
		std::transform(qstr.begin(), qstr.end(), qstr.begin(),
			[](char c) -> char {return c + 33;}); //qual value to actual str
		//returns 0 on success, -1 on fail. We should consider throwing if it fails.
		int r = bam_aux_update_str(this->r, "OQ", qstr.length(), qstr.c_str());
	}
	for(int i = 0; i < this->r->core.l_qseq; ++i){
		q[i] = (char)qual[i];
	}
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

void FastqFile::recalibrate(std::vector<int> qual){
	for(int i = 0; i < this->r->qual.l; ++i){
		this->r->qual.s[i] = (char)(qual[i]+33);
	}
}

int FastqFile::open_out(std::string filename){
	ofh = bgzf_open(filename.c_str(),"w"); //mode should be "wu" if uncompressed output is desired.
	return ofh == 0 ? -1 : 0;
}

int FastqFile::write(){
	std::string name = this->r->name.s ? std::string(this->r->name.s) : "";
	std::string seq = this->r->seq.s ? std::string(this->r->seq.s) : "";
	std::string comment = this->r->comment.s ? std::string(this->r->comment.s) : "";
	std::string qual = this->r->qual.s ? std::string(this->r->qual.s) : "";
	std::string s("@" + name + "\n" + seq + "\n+" + comment + "\n" + qual + "\n");
	return bgzf_write(ofh, s.c_str(), s.length());
}

// KmerSubsampler

//return the next hashed kmer
//once the file is finished iterating and there are no remaining kmers,
//return 0. That means you should check readseq.empty() if you get a 0!
uint64_t KmerSubsampler::next_kmer(){
	if(cur_kmer < kmers.size()){
		return kmers[cur_kmer++]; //return the current kmer and advance
	} else {
		readseq = file->next_str();
		if(readseq.empty()){
			return 0; //no more sequences
		} else {
			kmers = bloom::hash_seq(readseq, k); //get new vector of hashes
			cur_kmer = 0; //reset current kmer
			return this->next_kmer(); //try again
		}
	}
}

//return the next kmer that survives sampling
//once there are no more kmers return 0.
// check to see if readseq.empty() if you get a 0 result!
uint64_t KmerSubsampler::next(){
	uint64_t kmer = this->next_kmer();
	if(!readseq.empty()){
		if(d(rng)){ //sampled
			return kmer;
		}
		else{ //try again
			return this->next();
		}
	}
	return 0; // empty
}








}
