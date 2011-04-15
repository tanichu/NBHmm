//#include "../cpplapack_plus/cpplapack_plus.h"
#include "nbhmm.h"

int main(){
	
	scatter();

	
	int states = 3, dim=2;
	MLHmm H,Hest;	
	H.resize(states,dim);	
	Hest.resize(states,dim);
	
	H.read_TM("tm_.sample.dat");
	H.read_Mu("mu.sample.dat");
	H.read_diag_Sig("sigdiag.sample.dat");

	Hest.read_Mu("mu.sample.dat");
	Hest.read_diag_Sig("sigdiag.sample.dat");
	Hest.read_TM("tm.sample.dat");
	//cout << H.G[0].Mu<< endl;
	//cout << H.G[1].Mu<< endl;
	//cout << H.G[2].Mu<< endl;
	
	// Generation of Training Sample dataset
	vector<int> s  = GenerateStates(H,1000,0);
	//cout << dco(s) <<endl;
	
	dco(s).write("states.txt");
	dgematrix Y = GenerateObservations(H,s);
	
	//cout << Y ;
	Y.write("observation.txt");
	
	// Filtering by true hmm
	dgematrix F = ForwardFiltering(H,Y);
	F.write("ff.txt");
	dgematrix B = BackwardFiltering(H,Y);
	B.write("bf.txt");
	
	getchar();
	
	//Estimation of hmm
	int TRIAL = 10;
	dcovector likely_log(TRIAL);
	
	for (int i=0; i<TRIAL; i++) {
		
		//dcovector C;
		//ForwardBackwardFiltering(H,Y,F,B,C);
		
		
		//getchar();
		//cout << Hest.TM();
		//getchar();
		F = ForwardFiltering(Hest,Y);
		F.write("ffest.txt");
	/*	for (int t=0; t<F.m; t++) {
			double test=0;
			for (int j=0; j<F.n-1; j++) {
				test+=F(t,j);
			}
			cout << test << endl;
		}
		getchar();
	*/	
		B = BackwardFiltering(Hest,Y);
		B.write("bfest.txt");
		
		
	/*	for (int t=0; t<B.m; t++) {
			double test=0;
			for (int j=0; j<B.n-1; j++) {
				test+=B(t,j);
			}
			cout << test << endl;
		}
		
	*/	
		
		//cout << "Updating" <<endl;
		Hest.Update(Y,F,B);
		
		//cout << F(F.m-1,Hest.G.size())<<endl;
		//likely_log(i)=F(F.m-1,Hest.G.size());
		
		//cout << Hest.TM();
		
		//cout << sum_to_dro(Hest.TM_buffer);
		//cout << Hest.TM_buffer;
		
		getchar();
	}
	likely_log.write("lh.transition.txt");
	//dco(est).write("estimated.learned.txt");
	
	//cout << Hest.TM() << endl;
	
	


  return 0;

}
