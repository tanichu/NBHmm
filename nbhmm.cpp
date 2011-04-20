#include "nbhmm.h"




void NBHmm::read_TM(const char *filename){
	dgematrix mat;
	mat.read(filename);
	
	for (int i=0; i<M.size(); i++) {
		M[i].Mu = covec_read(mat,i);
	}
	
	TM();
	
}

// filename has each Mu as a drovector
// 
void NBHmm::read_Mu(const char *filename){
	dgematrix mat;
	mat.read(filename);
	
	for (int i=0; i<G.size(); i++) {
		G[i].Mu = t(rovec_read(mat,i));
	}
	
}



// filename has Sig's diag elements
// as drovectors 
void NBHmm::read_diag_Sig(const	char *filename){
	dgematrix mat;
	mat.read(filename);
	
	for (int i=0; i<G.size(); i++) {
		G[i].Sig = diag(t(rovec_read(mat,i)));	}
}

dgematrix NBHmm::TM(){
	dgematrix tm(M.size(),M.size());
	for (int i=0; i< M.size(); i++) {
		vec_set(tm,i,M[i].Mu);
	}
	
	TM_buffer = tm;
	return tm;
}

void NBHmm::resize(int n/*num_states*/,int d/*dim_output*/){
	//Gaussian
	G.resize(n);
	for (int i=0; i<n ; i++) {
		G[i].resize(d);
	}
	//Transition Matrix : Multinomial distribution 
	M.resize(n);
	for (int i=0; i<n ; i++) {
		M[i].resize(n);
	}
	pi.resize(n);
	all_set(pi,1.0/double(n));
	
};

// additional column is introduced
dgematrix ForwardFiltering(NBHmm H, dgematrix X){
	
	int num = H.M.size();
	
	double log_c_ = 0; 
	
	dgematrix F(X.m,num+1); // the last column is for baseline
	F.zero();
	
	dcovector current = H.TM() * H.pi;
	drovector x = rovec_read(X,0);	
	for (int i=0; i<num; i++) {
		current(i) = H.G[i].Probability(t(x));
	}
	double a = sum(current);
	current = (1.0/a)*current;
	log_c_ = log10(1.0/a);//similar to hamahata
	
	
	vec_set(F,0,t(direct_sum(current,log_c_)));
	
	
	for(int j=1;j<X.m;j++){
		current = H.TM_buffer*current;
		x = rovec_read(X,j);
		for (int i=0; i<current.l; i++) {
			current(i) = H.G[i].Probability(t(x))*current(i);
		}
		double a = sum(current);
		current = (1.0/a)*current;
		log_c_ = log10(1.0/a);//similar to hamahata
		vec_set(F,j,t(direct_sum(current,log_c_)));
	}
	
	return F;
}




vector<int> GenerateStates(NBHmm H, int length, int initial_state){
	
	vector<int> s; s.resize(length);
	s[0] = initial_state;
	for (int t=0; t<length-1; t++) {
		s[t+1] = H.M[s[t]].Sampler();
	}
	
	
	return s;
}


dgematrix GenerateObservations(NBHmm H, vector<int> s){
	
	dgematrix Y(s.size(),H.G[0].Mu.l);
	
	for (int i=0; i<Y.m; i++) {
		vec_set(Y,i,t(H.G[s[i]].Sampler()));
	}
	
	
	return Y;
	
}



void NBHmm::Update(dgematrix Y, vector<int> label){
	
	dgematrix TC = TransitionCount(label,M.size());
	int states = M.size();
	for (int i=0; i<states; i++) {
		//update Gaussians
		dgematrix Yi = ExtractMatrixByIndex(Y,label,i);
		G[i].UpdateMu(Yi);
		G[i].UpdateSig(Yi);
		//update Transition
		M[i].UpdateMu(covec_read(TC,i));
	}
	
}



vector<int> ForwardSampling(NBHmm H,dgematrix B, dgematrix X){
	
	
	vector<int> est;est.resize(B.m);
	
	dcovector current = subvector(t(rovec_read(B,0)),B.n-1);
	est[0] = MultiNominalSampler(current);
	
	dcovector f(H.M.size());
	
	for (int i=0; i<B.m-1; i++) {
		for (int k=0; k<B.n-1; k++) {
			dcovector x = t(rovec_read(X,i+1));
			f(k) = B(i+1, k)*H.M[est[i]].Probability(k)*H.G[k].Probability(x);
		}
		est[i+1] = MultiNominalSampler(f);
	}
	
	return est;
}



dgematrix BackwardFiltering(NBHmm H, dgematrix X){
	int num = H.M.size();
	
	double log_c_ = 0;
	
	dgematrix B(X.m, num+1);// the last column is for log baselie
	B.zero();
	// B(t,*) corresponds to m_{t+1,t}(*)
	dcovector current(num);
	
	all_set(current,1.0);
	vec_set(B,B.m-1,t(direct_sum(current,log_c_)));
	
	//B(T,*) is obviously {1} therefore we omit storing
	// the value to B(B.m-1,*)
	// So, B(B.m-1,*) contains beta(T-1)
	// index slides with 1 
	// This means F(0,*) and B(1,*) representing same time index 
	
	dgematrix tTM = t(H.TM());	
	drovector x;
	
	
	for (int j=X.m-1; j>0; j--) {//j is  time step
		
		x = rovec_read(X,j);
		for (int i=0; i<current.l; i++) {
			current(i) = H.G[i].Probability(t(x))*current(i);
		}
		
		current = tTM*current;		
		
		double a = sum(current);
		current = (1.0/a)*current;
		log_c_ = log10(1.0/a);//similar to hamahata
		
		
		vec_set(B,j-1,t(direct_sum(current,log_c_)));
	}
	
	// this scaling is not same as F's scaling parametere
	// http://unicorn.ike.tottori-u.ac.jp/murakami/doctor/node15.html
	// When forward-backward algorithm is used, 
	// we should take care about this difference.
	return B;
	
}



void NBShdpHmm::resize(int n/*num_max_states*/, int d/*dim_output*/){
	
	beta.resize(n);
	all_set(beta,1.0/double(n));
	hp_kappa = 0.9;
	hp_alpha = 0.1;
	hp_gamma = 0.1;
	
	//Gaussian
	G.resize(n);
	for (int i=0; i<n ; i++) {
		G[i].resize(d);
	}
	//Transition Matrix : Multinomial distribution 
	M.resize(n);
	for (int i=0; i<n ; i++) {
		M[i].resize(n);
	}
	
	
}



double NBShdpHmm::hp_rho(){
	return hp_kappa/(hp_alpha+hp_kappa);
}


void NBShdpHmm::Update_shdp(dgematrix Y,vector<int> label){
	int states = M.size();
	
	//Monitoring sampled Transition Count
	dgematrix TC = TransitionCount(label,M.size());	
	cout << "TC "<< TC << endl;
	cout << "augment variable m" << endl;
	
	//augment variable m
	dgematrix m(states,states);
	for (int j=0; j<states; j++) {
		for (int k=0; k<states; k++) {
			int n=0;
			for (int l=0; l<TC(k,j); l++) {
				//cout << "beta" << beta << endl;
				double p = (hp_alpha*beta(k)+hp_kappa*Kronecker_delta(j,k))/(n+hp_alpha*beta(k)+hp_kappa*Kronecker_delta(j,k));
				m(k,j) += BernoulliSampler(p);
				n++;
			}
		}
	}
	
	cout << "overide variable w " << endl;
	dcovector w(states);// overide variable
	for (int j=0; j<states; j++) {
		w(j) = BinominalSampler(m(j,j),hp_rho()/(hp_rho()+beta(j)*hp_rho()));
	}
	
	cout << "m_" << endl;
	dgematrix m_(states,states);
	for(int j=0;j<states;j++){
		for (int k=0; k<states; k++) {
			m_(k,j) = m(k,j) - w(j)*Kronecker_delta(j,k);
		}
	}
	
	//cout << "beta_prior" << endl;
	
	dcovector beta_prior(states);
	dcovector m_sum = sum_to_dco(m_);
	for (int i=0; i<states; i++) {
		beta_prior(i) = hp_gamma/states + m_sum(i);
	}
	
	//for lower bound of beta element (heuristics )
	for (int i=0; i < beta.l; i++) {
		if (beta(i)<DIR_MIN) {
			beta(i)=DIR_MIN;
			//cout << "dir lower cut!!"<< endl;
			//getchar();
		}
	}
		
	beta = DirichletSampler(beta_prior);	
	
	for (int i=0; i<states; i++) {
		//cout << i <<"-th module" << endl;
		//update Gaussians
		//cout << "ud gaussian "<<endl;
		dgematrix Yi = ExtractMatrixByIndex(Y,label,i);
		G[i].UpdateMu(Yi);
		G[i].UpdateSig(Yi);
		//update Transition
		//cout << "ud transition"  <<endl;
		M[i].hp_alpha = hp_alpha*beta;
		M[i].hp_alpha(i) += hp_kappa;
		M[i].UpdateMu(covec_read(TC,i));
	}
	
}


void MLHmm::Update_bw(dgematrix Y,dgematrix F, dgematrix B){
	//In Baum-Welch algorithm Forward and Backward filters
	// have to share the same FF'sscaling parameter 
	TM();
	
	int states = G.size();
	
	dcovector fC(F.m),bC(B.m);//buffer for scaling parameters
	
	//Buffering cumulative log forward scaling parameter 
	fC(F.m-1)=F(F.m-1,states);
	for(int k =F.m-1; k>0; k--){
		fC(k-1) = fC(k)+F(k-1,states);
	}
	
	//Buffering cumulative log backward scaling parameter
	bC(B.m-1)=B(B.m-1,states);
	for(int t=B.m-1; t>0; t--) {
		bC(t-1) = bC(t)+B(t-1,states);
	}
	
	
	//Preparing B_
	//b/pow(10.0,bC(t)) が真のb これに pow(10.0,fC(t))をかけると
	//Ffでスケーリングしたb'
	dgematrix B_(B.m,states);
	for(int t=0;t<B.m;t++){
		for (int i=0; i<states; i++) {
			B_(t,i)=B(t,i)*pow(10.0,fC(t)-bC(t));
		}
	}
	
	cout <<"before for loop" <<endl;
	dgematrix trans(3,3);
	
	vector<dgematrix> eta;
	eta.resize(Y.m);
	for (int t=0; t<Y.m;t++) {
		eta[t].resize(states,states);
	}
	for (int i=0; i<states; i++) {
		for (int j=0; j<states; j++) {
			dcovector upvec(states); upvec.zero();// = 0;
			double up = 0;
			double down = 0;
			for (int t=0; t<Y.m-1; t++) {
				
				eta[t](i,j)= F(t,i)*TM_buffer(j,i)*G[j].Probability(CPPL::t(rovec_read(Y,t+1)))*B_(t+1,j);
			}
		}		
	}
	
	//Calculating Pr(O|¥lambda)
	double Pr = 0;
	for (int i=0; i<states; i++) {
		for(int j=0; j<states; j++){
			Pr += eta[1](i,j);
		}
	}
	for (int t=0; t<Y.m; t++) {
		eta[t]= (1.0/Pr)*eta[t];
	}
	dgematrix gamma(Y.m,states);gamma.zero();
	for (int t=0; t<Y.m; t++) {
		for(int i=0;i<states;i++){
			for (int j=0; j<states; j++) {
				gamma(t,i) += eta[t](i,j);
			}
			
		}
	}
	
	dgematrix trans_(states,states);
	for(int i=0;i<states;i++){
		for (int j=0; j<states; j++) {
			double up=0;double down=0;
			for (int t=0; t<Y.m; t++) {
				up += eta[t](i,j);
				down += gamma(t,i);
			}
			trans_(j,i) = up/down;	
			M[i].Mu(j) = up/down;
		}
	}
	
	
	//	cout << "est transition" << endl;
	//	cout << trans_ << endl;
	//	cout << " transition " << (TM());
	//	cout << "sum transition " << sum_to_dro(TM());
	
	for (int i=0; i<states; i++) {
		cout << "state: "<< i << endl;
		cout << CPPL::t(G[i].Mu);
		cout << G[i].Sig;
		cout << M[i].Mu;
	}
	cout << "before " << endl;
	
	
	// Mu estimation
	for (int i=0; i<states; i++) {
		dcovector temp(Y.n);temp.zero();
		double probcount = 0;
		for(int t=0;t<Y.m;t++){
			probcount += gamma(t,i);
			temp = temp + gamma(t,i)*CPPL::t(rovec_read(Y,t));
		}
		G[i].Mu = (1.0/probcount)*temp;
	}
	
	// Sig estimation
	for (int i=0; i<states; i++) {
		dgematrix temp(Y.n,Y.n);temp.zero();
		double probcount = 0;
		for(int t=0;t<Y.m;t++){
			
			probcount += gamma(t,i);
			drovector tempro = rovec_read(Y,t);
			tempro = tempro - CPPL::t(G[i].Mu);
			temp = temp + gamma(t,i)*(CPPL::t(tempro)* tempro) ;
		}
		G[i].Sig = (1.0/(probcount+ G[i].hp_n))*(temp + G[i].hp_A);
	}
	
	
	// Update log show
	for (int i=0; i<states; i++) {
		cout << "state: "<< i << endl;
		cout << CPPL::t(G[i].Mu);
		cout << G[i].Sig;
	}
	
	cout << " transition " << endl << (TM());
	
	
}

vector<int> BackwardSampling(NBHmm H,dgematrix F){
	
	
	vector<int> est;est.resize(F.m);
	
	dcovector current = subvector(t(rovec_read(F,F.m-1)),F.n-1);
	est[F.m-1] = MultiNominalSampler(current);
	
	for (int i=F.m-1; i>0; i--) {
		for (int k=0; k<F.n-1; k++) {
			current(k) = F(i-1, k)*H.M[k].Probability(est[i]);
		}
		est[i-1] = MultiNominalSampler(current);
	}
	
	return est;
}

