//#include "../cpplapack_plus/cpplapack_plus.h"
#include "nbhmm.h"



int main(int argc, char *argv[]){
	
  srand(0);
	//scatter();	

	int states = 3, dim=2;
	NBHmm H;
	NBShdpHmm Hest;	
	H.resize(states,dim);	
	Hest.resize(2*states,dim);
	
	cout << Hest.beta.l << endl;
	cout << Hest.beta << endl;

	cout << "1"<<endl;
	H.read_TM("sample_dat/tm.sample.dat");
	H.read_Mu("sample_dat/mu.sample.dat");
	H.read_diag_Sig("sample_dat/sigdiag.sample.dat");

	
	cout << "2"<<endl;
	
	cout << "Start generating data set." << endl;
	
	// Generation of Training Sample dataset
	
	vector< vector<int> > ss;
	for (int i=0; i<3; i++) {
		ss.push_back(GenerateStates(H,100,0));	
	}
	dco(ss[0]).write("states_0.txt");
	dco(ss[1]).write("states_1.txt");
	dco(ss[2]).write("states_2.txt");
	
	vector< dgematrix > sY;
	
	sY.resize(3);
	sY[0] = GenerateObservations(H,ss[0]);	
	sY[1] = GenerateObservations(H,ss[1]);	
	sY[2] = GenerateObservations(H,ss[2]);	

	
	sY[0].write("observation_0.txt");
	sY[1].write("observation_1.txt");
	sY[2].write("observation_2.txt");

	
	
	//getchar();
	
	cout << "3"<<endl;
	
	
	// Filtering by true hmm
/*	dgematrix B= BackwardFiltering(H,Y);
	B.write("bf.txt");
	vector<int> est = ForwardSampling(H,B,Y);
	dco(est).write("estimated.txt");
	cout << "End" << endl;
*/
	
	cout << "Start estimation by using sticky HDP-HMM     ....press key!" << endl;
	getchar();
	
	//Estimation of hmm
	int TRIAL = 100;
	dcovector likely_log(TRIAL);

	cout << "4" << endl;

	
	vector< vector< int > > est;

	for (int i=0; i<TRIAL; i++) {
		est.resize(3);
		for (int j=0; j<3; j++) {
			dgematrix B = BackwardFiltering(Hest,sY[j]);
			est[j] = ForwardSampling(Hest,B,sY[j]);	
		}
		Hest.Update_shdp_multi(sY,est);
		
		//cout << B(0,Hest.G.size())<<endl;
		//likely_log(i)=B(0,Hest.G.size());
	}
	//likely_log.write("lh.transition.txt");
	//dco(est[0]).write("estimated.learned.txt");
	
	//cout << Hest.TM() << endl;
	
	
	for (int i=0; i<2*states; i++) {
		cout << "state: "<< i << endl;
		cout << t(Hest.G[i].Mu);
		cout << Hest.G[i].Sig;
		cout << Hest.M[i].Mu;
	}
	
	

	
  return 0;

}
