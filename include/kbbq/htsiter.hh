#ifndef KBBQ_HTSITER_H
#define KBBQ_HTSITER_H

#include <iterator>
#include <string>
#include <random>
#include <algorithm>
#include <htslib/hts.h>
#include <htslib/bgzf.h>
#include <htslib/sam.h>
#include <htslib/kseq.h>
#include <minion.hpp>
#include "readutils.hh"
#include "kseq.hh"

//This is defined starting in version 1.10
#ifndef HTS_VERSION
#define HTS_VERSION 0 
#endif

static_assert(HTS_VERSION >= 101000, "Your version of htslib is out of date. KBBQ requires version >= 10.1.");

//fwd declare
namespace readutils{
	class CReadData;
	std::string bam_seq_str(bam1_t* bamrecord);
}

//unsigned char* s = bam_get_seq(bamrecord);
namespace htsiter{

class HTSFile{
public:
	virtual ~HTSFile(){}
	virtual int next()=0;
	virtual std::string next_str()=0;
	virtual readutils::CReadData get()=0;
	virtual void recalibrate(std::vector<int> qual)=0;
	virtual int open_out(std::string filename)=0; //open an output file so it can be written to later.
	virtual int write()=0; //write the current read to the opened file.
};

class BamFile: public HTSFile{
public:
	samFile *sf;
	hts_idx_t *idx;
	sam_hdr_t *h;
	hts_itr_t *itr;
	bam1_t *r;
	samFile *of;
	bool use_oq;
	bool set_oq;
	BamFile(std::string filename, bool use_oq = false, bool set_oq = false):
		use_oq(use_oq), set_oq(set_oq){
		r = bam_init1();
		sf = sam_open(filename.c_str(), "r");
	   	idx = sam_index_load(sf, filename.c_str());
	    h = sam_hdr_read(sf);
	    itr = sam_itr_queryi(idx, HTS_IDX_START, 0, 0); //iterate over whole file;
	    of = NULL;
	}
	~BamFile(){
		if(r != NULL){bam_destroy1(r);}
		if(itr != NULL){sam_itr_destroy(itr);}
		if(sf != NULL){sam_close(sf);}
		if(of != NULL){sam_close(of);}
		if(h != NULL){sam_hdr_destroy(h);}
		if(idx != NULL){hts_idx_destroy(idx);}
	}

	// to use: while (ret = BamFile.next() >= 0){//do something with this->r}
	int next();
	// return next read sequence as a string. if there are no more, return the empty string.
	std::string next_str();
	//
	readutils::CReadData get();
	//
	void recalibrate(std::vector<int> qual);
	// TODO:: add a PG tag to the header
	int open_out(std::string filename);
	//
	int write();
	//
}; //end of BamFile class

class FastqFile: public HTSFile
{
public:
	BGZF* fh;
	kseq::kseq_t* r;
	BGZF* ofh;
	FastqFile(std::string filename){
		fh = bgzf_open(filename.c_str(),"r");
		r = kseq::kseq_init(fh);
	};
	~FastqFile(){
		if(r != NULL){kseq_destroy(r);}
		if(fh != NULL){bgzf_close(fh);}
		if(ofh != NULL){bgzf_close(ofh);}
	}
	int next();
	std::string next_str();
	readutils::CReadData get();
	void recalibrate(std::vector<int> qual);
	int open_out(std::string filename);
	int write();
};

class KmerSubsampler{
public:
	HTSFile* file;
	minion::Random rng;
	std::bernoulli_distribution d;
	std::string readseq = "";
	std::vector<uint64_t> kmers;
	size_t cur_kmer = 0;
	int k;
	KmerSubsampler(HTSFile* file): KmerSubsampler(file, 32){}
	KmerSubsampler(HTSFile* file, int k): KmerSubsampler(file, k, .15){}
	KmerSubsampler(HTSFile* file, int k, double alpha): KmerSubsampler(file, k, alpha,  minion::create_seed_seq().GenerateOne()){}
	KmerSubsampler(HTSFile* file, int k, double alpha, uint64_t seed): file(file), k(k), d(alpha) {rng.Seed(seed);}
	
	//return the next hashed kmer
	//once the file is finished iterating and there are no remaining kmers,
	//return 0. That means you should check readseq.empty() if you get a 0!
	//if !readseq.empty() that means a kmer hashed to 0 but there are further kmers.
	uint64_t next_kmer();

	//return the next hashed kmer that survives sampling
	//once there are no more kmers return 0.
	// check to see if readseq.empty() if you get a 0 result!
	//if !readseq.empty() that means a kmer hashed to 0 but there are further kmers.
	uint64_t next();
};

}

#endif
