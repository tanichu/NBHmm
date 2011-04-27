#include "../cpplapack_plus/cpplapack_plus.h"




//NBHmm  Forwardfiltering - Backward sampling のBayseHMM
class NBHmm  {
public:
	vector<NBGauss> G;//Gaussian distribution
	vector<NBMulti> M;// Transition multinomial distribution
	dcovector pi;// Initial distribution. (Not used???)
	void resize(int num_states, int dim_output);
	dgematrix TM_buffer; // used to buffering TM for reducing repetedly estimation of TM.
	dgematrix TM();// shows transition matrix
	
	
	void Update(dgematrix Y,vector<int> label);
	
	void read_Mu(const char *filename);
	void read_diag_Sig(const char *filename);
	
	void read_TM(const char *filename);// read TM from file
};


//For bayse HMM FF-BS
dgematrix ForwardFiltering(NBHmm H, dgematrix X);
vector<int> BackwardSampling(NBHmm H,dgematrix F);

//Generative process
vector<int> GenerateStates(NBHmm H, int length, int initial_state);
dgematrix GenerateObservations(NBHmm H, vector<int> s);


// designed for sticky HDP-HMM Backward filtering 
vector<int> ForwardSampling(NBHmm H,dgematrix B, dgematrix X);
dgematrix BackwardFiltering(NBHmm H, dgematrix X);



class MLHmm : public NBHmm{
public:
	//vector<NBGauss> G;//Gaussian distribution
	//vector<NBMulti> M;// Transition multinomial distribution
	
	//void resize(int num_states, int dim_output);
	//dgematrix TM_buffer; // used to buffering TM for reducing repetedly estimation of TM.
	//dgematrix TM();// shows transition matrix
	
	//void Update_bw(dgematrix Y,dgematrix F/*Forward message*/, dgematrix B/*Backword message*/);
	void Update_bw(dgematrix Y);
	//void read_Mu(const char *filename);
	//void read_diag_Sig(const char *filename);
	
	//void read_TM(const char *filename);// read TM from file
};




//forward-backward はそのまま行けるが，
// Updateについては変更が必要,また，Dirichlet 分布のハイパーパラメータが常時更新される．
class NBShdpHmm : public MLHmm{
public:
	dcovector beta;//global transition
	
	//	vector<NBGauss> G;//Gaussian distribution
	//	vector<NBMulti> M;// Transition multinomial distribution
	
	void resize(int num_max_states, int dim_output);
	//	dgematrix TM_buffer; // used to buffering TM for reducing repetedly estimation of TM.
	//	dgematrix TM();// shows transition matrix
	
	
	void Update_shdp(dgematrix Y,vector<int> label);
	
	//use	void read_Mu(const char *filename);
	//use	void read_diag_Sig(const char *filename);
	//use	void read_TM(const char *filename);// read TM from file
	
	
	//hyper parameter
	double hp_kappa;
	double hp_rho();
	double hp_gamma;
	double hp_alpha;
};

double loglikelihood(dgematrix F);


#define HMMclass NBShdpHmm 

