#include "covariateutils.hh"

namespace covariateutils{

	std::vector<long double> NormalPrior::normal_prior{};

	long double NormalPrior::get_normal_prior(size_t j){
		if(j >= normal_prior.size()){
			for(int i = normal_prior.size(); i < j+1; ++i){
				errno = 0;
				// long double prior_linspace = .9l * std::exp(-std::pow((long double)i,2.0l) * 2.0l)
				normal_prior.push_back(std::log(.9l * std::exp(-(std::pow(((long double)i/.5l),2.0l))/2.0l)));
				if(errno != 0){ //if an underflow happens just set the prior to smallest possible #
					normal_prior[i] = std::numeric_limits<long double>::lowest();
				}
			}
		}
		return normal_prior[j];
	}

	void CCovariate::increment(size_t idx, covariate_t value){
		(*this)[idx][0] += value[0];
		(*this)[idx][1] += value[1];
	}

	void CRGCovariate::consume_read(const readutils::CReadData& read){
		int rg = read.get_rg_int();
		if(this->size() <= rg){this->resize(rg+1);}
		// std::vector<bool> nse = read.not_skipped_errors();
		for(size_t i = 0; i < read.skips.size(); ++i){
			if(!read.skips[i]){
				this->increment(rg, std::array<unsigned long long,2>({read.errors[i], 1}));
			}
		}
	}

	//TODO: fix priors
	rgdq_t CRGCovariate::delta_q(meanq_t prior){
		rgdq_t dq(this->size());
		for(int i = 0; i < this->size(); ++i){ //i is rgs here
			int map_q = 0; //maximum a posteriori q
			long double best_posterior = std::numeric_limits<long double>::lowest();
			for(int possible = 0; possible < MAXQ+1; ++possible){
				int diff = std::abs(prior[i] - possible);
				long double prior_prob = NormalPrior::get_normal_prior(diff);
				long double p = recalibrateutils::q_to_p(possible);
				long double loglike = log_binom_pmf((*this)[i][0] + 1, (*this)[i][1] + 2, p);
				long double posterior = prior_prob + loglike;
				if(posterior > best_posterior){
					map_q = possible;
					best_posterior = posterior;
				}
			}
			dq[i] = map_q - prior[i];
		}
		return dq;
	}

	void CQCovariate::consume_read(const readutils::CReadData& read){
		int rg = read.get_rg_int();
		int q;
		if(this->size() <= rg){this->resize(rg+1);}
		for(size_t i = 0; i < read.skips.size(); ++i){
			if(!read.skips[i]){
				q = read.qual[i];
				if((*this)[rg].size() <= q){(*this)[rg].resize(q+1);}
				(*this)[rg].increment(q, std::array<unsigned long long, 2>({read.errors[i], 1}));
			}
		}
	}

	qscoredq_t CQCovariate::delta_q(prior1_t prior){
		qscoredq_t dq(this->size());
		for(int i = 0; i < this->size(); ++i){ //i is rgs here
			dq[i].resize((*this)[i].size());
			for(int j = 0; j < (*this)[i].size(); ++j){ //j is q's here
				int map_q = 0; //maximum a posteriori q
				long double best_posterior = std::numeric_limits<long double>::lowest();
				for(int possible = 0; possible < MAXQ+1; possible++){
					int diff = std::abs(prior[i] - possible);
					long double prior_prob = NormalPrior::get_normal_prior(diff);
					long double p = recalibrateutils::q_to_p(possible);
					long double loglike = log_binom_pmf((*this)[i][j][0] + 1, (*this)[i][j][1] + 2, p);
					long double posterior = prior_prob + loglike;
					if(posterior > best_posterior){
						map_q = possible;
						best_posterior = posterior;
					}
				}
				dq[i][j] = map_q - prior[i];
			}
		}
		return dq;
	}

	void CCycleCovariate::consume_read(const readutils::CReadData& read){
		int rg = read.get_rg_int();
		int q;
		int cycle;
		if(this->size() <= rg){this->resize(rg+1);}
		size_t readlen = read.skips.size();
		for(size_t i = 0; i < readlen; ++i){
			if(!read.skips[i]){
				q = read.qual[i];
				if((*this)[rg].size() <= q){(*this)[rg].resize(q+1);}
				if((*this)[rg][q][read.second].size() <= i){(*this)[rg][q][!!read.second].resize(i+1);}
				(*this)[rg][q][read.second].increment(i, std::array<unsigned long long, 2>({read.errors[i], 1}));
			}
		}
	}

	cycledq_t CCycleCovariate::delta_q(prior2_t prior){
		cycledq_t dq(this->size()); //rg -> q -> fwd(0)/rev(1) -> cycle -> values
		for(int i = 0; i < this->size(); ++i){ //i is rgs here
			dq[i].resize((*this)[i].size());
			for(int j = 0; j < (*this)[i].size(); ++j){ //j is q's here
				for(int k = 0; k < 2; ++k){ //fwd/rev
					dq[i][j][k].resize((*this)[i][j][k].size());
					for(int l = 0; l < (*this)[i][j][k].size(); ++l){ //cycle value
						int map_q = 0; //maximum a posteriori q
						long double best_posterior = std::numeric_limits<long double>::lowest();
						for(int possible = 0; possible < MAXQ+1; possible++){
							int diff = std::abs(prior[i][j] - possible);
							long double prior_prob = NormalPrior::get_normal_prior(diff);
							long double p = recalibrateutils::q_to_p(possible);
							long double loglike = log_binom_pmf((*this)[i][j][k][l][0] + 1, (*this)[i][j][k][l][1] + 2, p);
							long double posterior = prior_prob + loglike;
							if(posterior > best_posterior){
								map_q = possible;
								best_posterior = posterior;
							}
						}
						dq[i][j][k][l] = map_q - prior[i][j];
					}
				}
			}
		}
		return dq;
	}

	void CDinucCovariate::consume_read(const readutils::CReadData& read, int minscore){
		int rg = read.get_rg_int();
		if(this->size() <= rg){this->resize(rg+1);}
		int q;
		for(size_t i = 1; i < read.seq.length(); i++){
			q = read.qual[i];
			//not skipped guarantees read.seq[i] != 'N' and q >= minscore
			if(!read.skips[i] && seq_nt4_table[read.seq[i-1]] < 4){
				if((*this)[rg].size() <= q){(*this)[rg].resize(q+1);}
				if((*this)[rg][q].size() < 16){(*this)[rg][q].resize(16);}
				(*this)[rg][q].increment(dinuc_to_int(read.seq[i-1], read.seq[i]),std::array<unsigned long long, 2>({read.errors[i], 1}));
			}
		}
		// seq_nt16_table[256]: char -> 4 bit encoded (1/2/4/8)
		// seq_nt16_str[]: 4 bit -> char
		// seq_nt16_int[]: 4 bit -> 2 bits (0/1/2/3)
	}

	dinucdq_t CDinucCovariate::delta_q(prior2_t prior){
		dinucdq_t dq(this->size()); //rg -> q -> fwd(0)/rev(1) -> cycle -> values
		for(int i = 0; i < this->size(); ++i){ //i is rgs here
			dq[i].resize((*this)[i].size());
			for(int j = 0; j < (*this)[i].size(); ++j){ //j is q's here
				dq[i][j].resize((*this)[i][j].size());
				for(int k = 0; k < (*this)[i][j].size(); ++k){ //k is dinuc
					int map_q = 0; //maximum a posteriori q
					long double best_posterior = std::numeric_limits<long double>::lowest();
					for(int possible = 0; possible < MAXQ+1; possible++){
						int diff = std::abs(prior[i][j] - possible);
						long double prior_prob = NormalPrior::get_normal_prior(diff);
						long double p = recalibrateutils::q_to_p(possible);
						long double loglike = log_binom_pmf((*this)[i][j][k][0] + 1, (*this)[i][j][k][1] + 2, p);
						long double posterior = prior_prob + loglike;
						if(posterior > best_posterior){
							map_q = possible;
							best_posterior = posterior;
						}
					}
					dq[i][j][k] = map_q - prior[i][j];
				}
			}
		}
		return dq;
	}

	void CCovariateData::consume_read(readutils::CReadData& read, int minscore){
		for(int i = 0; i < read.seq.length(); ++i){
			read.skips[i] = (read.skips[i] || seq_nt4_table[read.seq[i]] >= 4); // || read.qual[i] < minscore);
		}
		rgcov.consume_read(read);
		qcov.consume_read(read);
		cycov.consume_read(read);
		dicov.consume_read(read, minscore);
	}

	dq_t CCovariateData::get_dqs(){
		dq_t dq;
		std::vector<long double> expected_errors(this->qcov.size(),0);
		meanq_t meanq(this->qcov.size(),0);
		for(int rg = 0; rg < this->qcov.size(); ++rg){
			for(int q = 0; q < this->qcov[rg].size(); ++q){
				expected_errors[rg] += (recalibrateutils::q_to_p(q) * this->qcov[rg][q][1]);
			}
			meanq[rg] = recalibrateutils::p_to_q(expected_errors[rg] / this->rgcov[rg][1]);
		}
		dq.meanq = meanq;
		dq.rgdq = this->rgcov.delta_q(meanq);
		prior1_t rgprior(this->qcov.size());
		for(int rg = 0; rg < this->qcov.size(); ++rg){
			rgprior[rg] = meanq[rg] + dq.rgdq[rg];
		}
		dq.qscoredq = this->qcov.delta_q(rgprior);
		prior2_t qprior(this->qcov.size());
		for(int rg = 0; rg < this->qcov.size(); ++rg){
			for(int q = 0; q < this->qcov[rg].size(); ++q){
				qprior[rg].push_back(rgprior[rg] + dq.qscoredq[rg][q]);
			}
		}
		dq.cycledq = this->cycov.delta_q(qprior);
		dq.dinucdq = this->dicov.delta_q(qprior);
		return dq;
	}



}