//#include "../cpplapack_plus/cpplapack_plus.h"
#include "nbhmm.h"

int main(){
	
	scatter();

	
	int states = 3, dim=2;
	MLHmm H,Hest;	
	H.resize(states,dim);	
	Hest.resize(states,dim);
	
	H.read_TM("sample_dat/tm.sample.dat");
	H.read_Mu("sample_dat/mu.sample.dat");
	H.read_diag_Sig("sample_dat/sigdiag.sample.dat");

	for(int i=0;i<Hest.G.size();i++){
		Hest.G[i].Mu = MultiGaussSampler(Hest.G[i].hp_m,Hest.G[i].hp_S);
		Hest.G[i].Sig.identity();
	}
	
	//Hest.read_Mu("sample_dat/mu_6.sample.dat");
	//Hest.read_TM("sample_dat/tm_6.sample.dat");
	
	printf("H Hest prepared\n");
	
	// Generation of Training Sample dataset
	vector<int> s  = GenerateStates(H,1000,0);
	//cout << dco(s) <<endl;
	dco(s).write("states.txt");
	
	
	dgematrix Y = GenerateObservations(H,s);
	//cout << Y ;
	Y.write("observation.txt");
	
	//Estimation of hmm
	int TRIAL = 100;
	dcovector likely_log(TRIAL);
	
	dgematrix F;
	double lk_=0;
	double lk=0;
	for (int i=0; i<TRIAL; i++) {
		Hest.Update_bw(Y);
		F = ForwardFiltering(Hest,Y);
		F.write("ffest.txt");
		
		cout << lk << endl;
	
		lk = sum_to_dro(F)(F.n-1);
		
		likely_log(i) = lk;
		
		getchar();
		
	}
	likely_log.write("lh.transition.txt");
	//dco(est).write("estimated.learned.txt");
	
	//cout << Hest.TM() << endl;
	
	


  return 0;

}
