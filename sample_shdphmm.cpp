//#include "../cpplapack_plus/cpplapack_plus.h"
#include "nbhmm.h"



int main(int argc, char *argv[]){
	
  srand(atoi(argv[1]));
	//scatter();	

	int states = 3, dim=2;
	NBHmm H;
	NBShdpHmm Hest;	
	H.resize(states,dim);	
	Hest.resize(2*states,dim);
	
	cout << Hest.beta.l << endl;
	cout << Hest.beta << endl;
	getchar();
	cout << "1"<<endl;
	H.read_TM("tm.sample.dat");
	H.read_Mu("mu.sample.dat");
	H.read_diag_Sig("sigdiag.sample.dat");

	
	cout << "2"<<endl;

	
	// Generation of Training Sample dataset
	vector<int> s  = GenerateStates(H,1000,0);	
	dco(s).write("states.txt");
	dgematrix Y = GenerateObservations(H,s);	
	Y.write("observation.txt");

	cout << "3"<<endl;
	
	
	// Filtering by true hmm
	dgematrix B= BackwardFiltering(H,Y);
	B.write("bf.txt");
	vector<int> est = ForwardSampling(H,B,Y);
	dco(est).write("estimated.txt");


	//Estimation of hmm
	int TRIAL = 100;
	dcovector likely_log(TRIAL);

	cout << "4" << endl;
	
	for (int i=0; i<TRIAL; i++) {
	//	cout << "for loop: " << i << endl;
	//	cout << "BF" <<endl;
		B = BackwardFiltering(Hest,Y);
	//	cout << "FS" <<endl;		
		est = ForwardSampling(Hest,B,Y);	
	//	cout << "Updating" <<endl;
		Hest.Update(Y,est);
		
		cout << B(0,Hest.G.size())<<endl;
		likely_log(i)=B(0,Hest.G.size());
	}
	likely_log.write("lh.transition.txt");
	dco(est).write("estimated.learned.txt");
	
	//cout << Hest.TM() << endl;
	
	
	for (int i=0; i<states; i++) {
		cout << "state: "<< i << endl;
		cout << t(Hest.G[i].Mu);
		cout << Hest.G[i].Sig;
		cout << Hest.M[i].Mu;
	}
	
	
	getchar();

	
  return 0;

}
