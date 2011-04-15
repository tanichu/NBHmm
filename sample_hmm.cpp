#include "nbhmm.h"

int main(int argc, char *argv[]){
	srand(atoi(argv[1]));
	//scatter();

	int big_sates = 6;
	int states = 3, dim=2;
	NBHmm H,Hest;	
	H.resize(states,dim);	
	Hest.resize(big_sates/*states*/,dim);
	
	H.read_TM("tm.sample.dat");
	H.read_Mu("mu.sample.dat");
	H.read_diag_Sig("sigdiag.sample.dat");
	
	cout << H.G[0].Mu<< endl;
	cout << H.G[1].Mu<< endl;
	cout << H.G[2].Mu<< endl;
	
	// Generation of Training Sample dataset
	vector<int> s  = GenerateStates(H,1000,0);
	cout << dco(s) <<endl;
	
	dco(s).write("states.txt");
	dgematrix Y = GenerateObservations(H,s);
	
	cout << Y ;
	Y.write("observation.txt");
	
	// Filtering by true hmm
	dgematrix F = ForwardFiltering(H,Y);
	F.write("ff.txt");
	vector<int> est = BackwardSampling(H,F);
	dco(est).write("estimated.txt");

	
	//Estimation of hmm
	int TRIAL = 100;
	dcovector likely_log(TRIAL);
	
	for (int i=0; i<TRIAL; i++) {
		//cout << "FF" <<endl;
		F = ForwardFiltering(Hest,Y);
		//cout << "BS" <<endl;		
		est = BackwardSampling(Hest,F);	
		//cout << "Updating" <<endl;
		Hest.Update(Y,est);
		
		cout << sum_to_dro(F)(F.n-1)<<endl;
		likely_log(i)=sum_to_dro(F)(F.n-1);
		
		cout << TransitionCount(est,big_sates)<<endl;
		//cout << Hest.TM() <<endl;
		//getchar();
	}
	likely_log.write("lh.transition.txt");
	dco(est).write("estimated.learned.txt");
		
	
	for (int i=0; i<states; i++) {
		cout << "state: "<< i << endl;
		cout << t(Hest.G[i].Mu);
		cout << Hest.G[i].Sig;
		cout << Hest.M[i].Mu;
	}
	
	getchar();
	
	/*
	dgematrix M("M.dat");
	vector<int>  label;label.resize(6);
	label[0]=1; 
	label[1]=2; 
	label[2]=3; 
	label[3]=1; 
	label[4]=2; 
	label[5]=3; 

	cout << ExtractMatrixByIndex(M,label,5) << endl;

	NBGauss G;
	G.resize(2);
	
	
	dgematrix X(10,2);
	
	for (int i=0; i <10; i++) {
		vec_set(X,i,t(G.Sampler()));
	}

	cout << X;
	
	
	cout << G.Mu << endl;
	cout << G.Sig << endl;
	
	getchar();
	
	NBGauss GL;GL.resize(2);
	GL.Mu(0) =10;
	//GL.Sig(0,0) = 10;
	
	cout << GL.Mu << GL.Sig << endl;
	
	getchar();*/
/*	
	for (int i=0; i<10; i++){
		GL.UpdateMu(X);
		cout << GL.Mu;
		GL.UpdateSig(X);
		cout << GL.Sig;
		getchar();
	}
*/
	/*
	NBMulti T; T.resize(4);
	dcovector v(4);v.zero();v(0) = 3;
	cout << v << endl;
	
	for (int i=0; i < 10; i++) {
		T.UpdateMu(v);
		cout << T.Mu << endl;
		cout << T.Sampler() << endl;
		getchar();
	}

	*/
	
/*	dcovector temp(2);temp.zero();
	for (double t=0;t<10; t+=0.1) {
		temp(0)=t;
		cout << G.Probability(temp) << endl;
	}
*/	
	
  return 0;

}
