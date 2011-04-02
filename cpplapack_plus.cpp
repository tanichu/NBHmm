#include "cpplapack_plus.h"
#include <algorithm>

using namespace std;
std::string itos(int n)
{
	std::stringstream str_stream;
	str_stream << n;
	return str_stream.str();
}

std::string ftos(double n)
{
	std::stringstream str_stream;
	str_stream << n;
	return str_stream.str();
}

string ScanTillRet(string str){
	return str.substr(0,str.find("\n"));
}
string ScanTillSpace(string str){
	return str.substr(0,str.find(" "));
}
string ScanTillTab(string str){
	return str.substr(0,str.find("\t"));
}
string ScanTillComma(string str){
	return str.substr(0,str.find(","));
}


string TrimLeft( string arg ) {
	int n = int(arg.find_first_not_of("\n"));
	int t = int(arg.find_first_not_of("\t"));
	int s = int(arg.find_first_not_of(" "));

	for(;(abs(n-s)+abs(n-t)+abs(t-s))!=0;){
	
		arg = arg.substr(min(n,min(t,s))+1,arg.length());
	n = int(arg.find_first_not_of("\n"));
	t = int(arg.find_first_not_of("\t"));
	s = int(arg.find_first_not_of(" "));
	}

	return arg;
}

string TrimRight( string arg ) {

	int n = int(arg.find_last_not_of("\n"));
	int t = int(arg.find_last_not_of("\t"));
	int s = int(arg.find_last_not_of(" "));
	
	while((abs(n-s)+abs(n-t)+abs(t-s))!=0){
	
		arg = arg.substr(0,min(n,min(t,s))+1);
		
		n = int(arg.find_last_not_of("\n"));
		t = int(arg.find_last_not_of("\t"));
		s = int(arg.find_last_not_of(" "));
	}


    return arg;
}



string Trim( string arg){
	return TrimRight(TrimLeft(arg));
}

string DelTillRet(string str){
	return str.substr(str.find_first_of("\n")+1,str.length());
}
string DelTillSpace(string str){
	return str.substr(str.find_first_of(" ")+1,str.length());
}
string DelTillTab(string str){
	return str.substr(str.find_first_of("\t")+1,str.length());
}

string DelTillComma(string str){
	return str.substr(str.find_first_of(",")+1,str.length());
}


vector<string> split(const string& str, const string& delimiter) {
    // delimiter(2 文字以上も可) を空白に置換
    std::string item(str);    
    for(unsigned pos = item.find(delimiter); pos != string::npos; pos = item.find(delimiter, pos)) {
        item.replace(pos, delimiter.size(), " ");
    }
    // 分解
    stringstream buf(item);

    // 読み取り
    std::vector<string> result;
    while(buf >> item) {
        result.push_back(item);
    }

    return result;
}


//csvを読み取る．
vector< vector<string> > csv_reader(const char *filename){
	char ch[10000];
	string str;
	ifstream fin(filename);
	vector< vector<string> > csv;

	while(fin.getline(ch,10000)){//データの読み込み
		str = ch;
		vector<string> strs;
		strs = split(str,",");
		csv.push_back(strs);
	}
	fin.close();

	return csv;
}

//ssvを読み取る．space separated value
vector< vector<string> > ssv_reader(const char *filename){
	char ch[10000];
	string str;
	ifstream fin(filename);
	vector< vector<string> > csv;

	while(fin.getline(ch,10000)){//データの読み込み
		str = ch;
		vector<string> strs;
		strs = split(str," ");
		csv.push_back(strs);
	}
	fin.close();

	return csv;
}	

	
///元suplapack

void scatter(){
      srand( (unsigned)time( NULL ) );
      setall(rand(),rand());
}


double gauss_rand(){
        
        double x1 = double(rand()+1.0)/(RAND_MAX +1.0);
        double x2 = double(rand()+1.0)/(RAND_MAX +1.0);

        double y = sqrt(-2*log(x1))*cos(2.0*M_PI*x2) ;
        return y;
}

double gauss_rand(double ave,double sig){
        return ave + sig*gauss_rand();
}
dcovector gauss_rand(dcovector ave,double sig){

        dcovector out(ave);
        int i;

        for(i=0;i<ave.l;i++){
                out(i) =gauss_rand(ave(i),sig) ;
        }
        return out;
}

dgematrix gauss_rand(dgematrix ave,double sig){

        dgematrix out(ave);
        int i,j;

        for(i=0;i<ave.m;i++){
        for(j=0;j<ave.n;j++){
                out(i,j) =gauss_rand(ave(i,j),sig) ;
		}}
        return out;
}


dcovector gauss_rand(dcovector ave,dcovector sig){

        dcovector out(ave);
        int i;

		for(i=0;i<ave.l;i++){
                out(i) =gauss_rand(ave(i),sig(i)) ;
        }
        return out;
}


double LGF(dcovector x,dcovector m,dcovector s){

        double log_p = 0;
        for(int i=0;i<x.l;i++){
			log_p += -0.5*(x(i)-m(i))*(x(i)-m(i))/(s(i)*s(i)) - log(sqrt(2*M_PI)*s(i)) ;
        }
        return log_p;
}

dcovector ng_calculation(dcovector *m,dcovector *s,int N,dcovector x){
//normalized gaussian
                dcovector lg(N),lg_(N),ng(N);

                double maxlg=0;
                double sumg_=0;

                for(int j=0;j<N;j++){
                       lg(j) =LGF(x,m[j],s[j]);
                       if((maxlg<lg(j))||(j==0)){maxlg=lg(j);}
                }
                for(int j=0;j<N;j++){
                        lg_(j) = lg(j)-maxlg;
                }
                sumg_=0;
                for(int j=0;j<N;j++){
                        sumg_ += exp(lg_(j));
                }
                for(int j=0;j<N;j++){
                        ng(j) = exp(lg_(j))/sumg_;
                }

                return ng;
}

dcovector rbf_calculation(dcovector *m,dcovector *s,int N,dcovector x){
//radial basis function gaussian
                dcovector lg(N),rbf(N);
               //argmaxlg=0;
                for(int j=0;j<N;j++){
                       lg(j) =LGF(x,m[j],s[j]);
                       rbf(j) = exp(lg(j));
                 }
                return rbf;
}

void  each_abs(dcovector &v){
        int i;
        for(i=0;i<v.l;i++){
                v(i) =fabs(v(i));
        }
}

void  each_square(dcovector &v){
        int i;
        for(i=0;i<v.l;i++){
                v(i) =v(i)*v(i);
        }
}

double frobnorm(dgematrix A){
    double norm = 0;     
	int i,j;
    for(i=0;i<A.m;i++){
        for(j=0;j<A.n;j++){
			norm+= A(i,j)*A(i,j);
		}
	}
	return sqrt(norm);
	
}







void restrict(dgematrix& A,double min,double max){
        int i,j;
        for(i=0;i<A.m;i++){
        for(j=0;j<A.n;j++){
                if(A(i,j)<min){A(i,j)=min;}
                if(A(i,j)>max){A(i,j)=max;}
        }}
}

void restrict(dcovector &v, double min,double max){

        int i;
        for(i=0;i<v.l;i++){
                if(v(i)<min){v(i)=min;}
                if(v(i)>max){v(i)=max;}
        }
}

void restrict(drovector &v, double min,double max){

        int i;
        for(i=0;i<v.l;i++){
                if(v(i)<min){v(i)=min;}
                if(v(i)>max){v(i)=max;}
        }
}

double trace(dgematrix A){

        int i;
        double trace=0;
        if(A.n==A.m){
                for(i=0;i<A.m;i++){
                       trace = A(i,i) + trace;
                }
        }
        else{
               cout << "正方行列でないのでトレースをとれません"<<endl;
        }
        return trace;
}

void read_multi_vector(dcovector *x,int num, int dim, char *filename){//一行10000バイトまでです．
	ifstream fin(filename);

	char ch[10000];
	string str;

	for(int i=0;i<num;i++){
		x[i].resize(dim);
		fin.getline(ch,10000);
		str = ch;
		for(int j=0;j<dim;j++){
			string xs = ScanTillSpace(str);
			str = DelTillSpace(str);
			str = TrimLeft(str);
			x[i](j) = atof(xs.c_str());
		}
	}

	fin.close();
}
void read_multi_vector(drovector *x,int num, int dim, char *filename){//一行10000バイトまでです．
	ifstream fin(filename);

	char ch[10000];
	string str;

	for(int i=0;i<num;i++){
		x[i].resize(dim);
		fin.getline(ch,10000);
		str = ch;
		for(int j=0;j<dim;j++){
			string xs = ScanTillSpace(str);
			str = DelTillSpace(str);
			str = TrimLeft(str);
			x[i](j) = atof(xs.c_str());
		}
	}

	fin.close();
}

dgematrix product(dcovector c,drovector r){
	dgematrix A(c.l,r.l),C(c.l,1),R(1,r.l);
	
	for(int i=0;i<c.l;i++){
		C(i,0)=c(i);
	}
	for(int j=0;j<r.l;j++){
		R(0,j)=r(j);
	}
	A = C*R;

	return A;
}



void vec_set(dgematrix &A, int k , dcovector x){
	if(k>=A.n){cout << "k is too learge for matrix A" <<endl;}
	else if(A.m!=x.l){cout << "vec dimension is not match to matrix's dimension."<<endl;}
	else{
		for(int i=0;i<A.m;i++){
			A(i,k)=x(i);
		}
	}
}
void vec_set(dgematrix &A, int k , drovector x){
	if(k>=A.m){cout << "k is too learge for matrix A" <<endl;}
	else if(A.n!=x.l){cout << "vec dimension is not match to matrix's dimension."<<endl;}
	else{
		for(int i=0;i<A.n;i++){
			A(k,i)=x(i);
		}
	}
}

void read_multi_matrix(char *filename,dgematrix *A,int n){
//read matrix files from "./filename/??.dat"
// ?? means indices from 0.
	string str,dir;
	str = filename;
	dir = string("./")+str;
	dir += "/";

	string head;
	head = dir;

	for(int i=0;i<n;i++){
		string file;
		file = head + itos(i)+string(".dat");
		A[i].read(file.c_str());
	}
}

int ProbSelect(dcovector v){
        dcovector p = v;

        double sum=0;
        for(int i=0;i<v.l;i++){
                sum+=v(i);
        }

        p *= (1.0/sum);

        double t = double(rand())/(double(RAND_MAX) + 1.0);
        int k =0;
        for(int i=0;i<v.l;i++){
                if(t>=0){
                        t=t-p(i);
                        if(t<0){k=i;}
                }
        }

        return k;
}

double max_value(dcovector v){
       double max =  v(0);
       for(int i=0;i<v.l;i++){if(v(i)>max){max = v(i);}}
       return max;
}

int argmax(dcovector v){

       int i;
       double max =  v(0);
       int argmax = 0;
       for(i=0;i<v.l;i++){
               if(v(i)>max){max = v(i);argmax=i;}
       }

       return argmax;

}

int argmin(dcovector v){
        int i;
        double min = v(0);
        int argmin = 0;
		for(i=0;i<v.l;i++){
                if(v(i)<min){min = v(i); argmin=i;}
        }

        return argmin;
}

double min_value(dcovector v){
        int i;
        double min = v(0);
        int argmin = 0;
        for(i=0;i<v.l;i++){
                if(v(i)<min){min = v(i); argmin=i;}
        }
        return min;
}

dcovector argmax_vector(dcovector v){
	dcovector w(v);
	w.zero();

	w(argmax(v)) = 1.0;
	return w;

}
dcovector argmin_vector(dcovector v){
	dcovector w(v);
	w.zero();

	w(argmin(v)) = 1.0;
	return w;
}



dcovector e_greedy(dcovector a,double greedy_ratio){
	double other_ratio=(1.0-greedy_ratio)/(double(a.l)-1.0);
	dcovector p(a.l);
	for(int i=0;i<a.l;i++){
		p(i)=other_ratio;
	}

	p(argmax(a)) = greedy_ratio;

	return p;

}

dcovector normalize(dcovector x){
	double sum = 0;
	for(int i=0;i<x.l;i++){
		sum += x(i);
	}

	return (1.0/double(sum))*x;
}

dcovector each_product(dcovector a,dcovector b){
	dcovector c =a;

	for(int i=0;i<a.l;i++){
		c(i) = a(i)*b(i);
	}

	return c;
}

double det(dgematrix q){

	dgematrix V,D;
	dcovector s;

	if(0!=q.dgesvd(s,V,D)){cout << "det error!"<< endl << "q=" <<q << endl;/*getch()*/;}

	double d = 1.0;
	for(int i=0;i<s.l;i++){
		d *= s(i);
	}

	return d;
}

double gaussian(dcovector x, dcovector m,dgematrix q){
	dcovector y  = x-m;
	drovector yt = t(y);
	dgematrix qi = i(q);
	return exp(-0.5*yt*qi*y)/sqrt(pow(2.0*M_PI,(double)x.l)*det(q));
}

double LGF(dcovector x, dcovector m,dgematrix q){
	dcovector y  = x-m;
	drovector yt = t(y);
	dgematrix qi = i(q);

	return -0.5*yt*qi*y-log(sqrt(pow(2.0*M_PI,(double)x.l)*det(q)));
}

double gaussian(dcovector x, dcovector m,dgematrix iq,double detq){
	dcovector y  = x-m;
	drovector yt = t(y);
	return exp(-0.5*yt*iq*y)/sqrt(pow(2.0*M_PI,(double)x.l)*detq);
}
double LGF(dcovector x, dcovector m,dgematrix iq,double detq){
	dcovector y  = x-m;
	drovector yt = t(y);
	return -0.5*yt*iq*y-log(sqrt(pow(2.0*M_PI,(double)x.l)*detq));
}

string csv(dcovector x){
	string out("");
	string comma(",");
	for(int i=0;i<x.l;i++){
		out += comma + ftos(x(i));
	}
	return out;
}


string csv(drovector x){
	string out("");
	string comma(",");
	for(int i=0;i<x.l;i++){
		out += comma + ftos(x(i));
	}
	return out;
}

string csv(dgematrix x){
	string out("");
	string comma(",");
	for(int i=0;i<x.m;i++){
	for(int j=0;j<x.n;j++){
		out += comma + ftos(x(i,j));
	}}
	return out;
}

void self_add(dgematrix &A,dgematrix B){
	for(int i=0;i<A.m;i++){
		for(int j=0; j<A.n;j++){
			A(i,j) += B(i,j);
		}
	}
}
void self_add(dcovector &x,dcovector y){
	for(int i=0;i<x.l;i++){
		x(i) += y(i);
	}
}
void self_set(dgematrix &A,dgematrix B){
	for(int i=0;i<A.m;i++){
		for(int j=0; j<A.n;j++){
			A(i,j) = B(i,j);
		}
	}
}
void self_set(dcovector &x,dcovector y){
	for(int i=0;i<x.l;i++){
		x(i) = y(i);
	}
}

void all_set(dgematrix &A,double value){
	for(int i=0;i<A.m;i++){
		for(int j=0; j<A.n;j++){
			A(i,j) = value;
		}
	}
}

void all_set(dcovector &v,double value){
	for(int i=0;i<v.l;i++){
		v(i) = value;
	}	
}


dgematrix submatrix(dgematrix A, int L, int R, int T, int B){
	dgematrix S(B-T+1,R-L+1);

	for(int i=0;i<B-T+1;i++){
		for(int j=0;j<R-L+1;j++){
			S(i,j)=A(T+i,L+j);		
		}
	}
	return S;
}

dgematrix lean_on (dgematrix A, dgematrix B){
	dgematrix C(A.m,A.n+B.n);

	for(int i=0;i<A.m;i++){
		for(int j=0;j<A.n;j++){
			C(i,j)=A(i,j);
		}
	}
	for(int i=0;i<B.m;i++){
		for(int j=0;j<B.n;j++){
			C(i,A.n+j)=B(i,j);
		}
	}
	return C;

}
dgematrix  put_on (dgematrix A, dgematrix B){
	dgematrix C(A.m+B.m,A.n);

	for(int i=0;i<A.m;i++){
		for(int j=0;j<A.n;j++){
			C(i,j)=A(i,j);
		}
	}
	for(int i=0;i<B.m;i++){
		for(int j=0;j<B.n;j++){
			C(A.m+i,j)=B(i,j);
		}
	}

	return C;
}

dcovector covec_read(dgematrix A, int k){
  dcovector v(A.m);
  for(int i=0; i< A.m;i++){
    v(i)= A(i,k); 
  }
  return v;
}

drovector rovec_read(dgematrix A, int k){
  drovector v(A.n);
  for(int i=0;i<A.n;i++){
    v(i) = A(k,i);
  }
  return v;
}



/// Hamahata procedres




//--------------------------------------------------------
//-------------------------関数----------------------------
//--------------------------------------------------------
//行列の指定した列に列ベクトルを挿入
dgematrix Insert_mat(dgematrix mat,int col,dcovector vec)
{
	for(int i=0;i<(vec.l);i++)
		mat(i,col) = vec(i);
    
	return(mat);
}

//行列の中の指定した列ベクトルを抜き出す
dcovector Extract_vec(dgematrix mat,int col)
{
	dcovector x(mat.m);
	
	for(int i=0;i<mat.m;i++)
		x(i) = mat(i,col);
    
	return(x);
}

//ベクトルのノルムを計算する
double Cal_norm(dcovector vec)
{
	double ans=0;
	for(int h=0;h<vec.l;h++)
		ans += pow(vec(h),2.0);
    
	return(pow(ans,0.5));
}

//行列にベクトルを挿入して足し合わせる
dgematrix InsertAdd_VecToMat(dgematrix mat,int col,dcovector vec)
{
	for(int i=0;i<(vec.l);i++)
		mat(i,col) += vec(i);
	
	return(mat);
}


//正則行列かどうかのチェック（正則なら0を返す）
//正則行列かどうかのチェック（正則なら0を返す）
int check_Regularization(dgematrix mat)
{
	std::vector<double> wr,wi;
	mat.dgeev(wr,wi);
	for(int i=0;i<mat.m;i++)
		if(wr[i] < 0)
			return(1);
	
	return(0);
}


//多次元正規分布の尤度を計算
double Cal_MultiNormLikely(dcovector x,dcovector mu,dgematrix sig)
{
	dcovector y = x - mu;
	double tmp1,tmp2;
	double ans;
	
	tmp1 = pow(pow(2.0*PI,double(y.l)),0.5) * pow(det(sig),0.5);
	
	tmp2 = -0.5 * CPPL::t(y) * CPPL::i(sig) * y;
	
	ans = exp(tmp2) / tmp1;
	if(ans == 0)
		ans = 1.0 / pow(10.0,100.0);
	
	return(ans);
}

//コレスキー分解
dgematrix cholesky(dgematrix mat)
{
	dgematrix L = mat;
	double tmp;
	L.zero();
	
	for(int i=0;i<mat.m;i++)
    {
		//i>jの計算
		for(int j=0;j<i;j++)
        {
			tmp = mat(i,j);
			for(int k=0;k<j;k++)
				tmp -= L(i,k) * L(j,k);
			
			L(i,j) = tmp / L(j,j);
        }
		//i=jの計算
		tmp = mat(i,i);
		for(int k=0;k<i;k++)
			tmp -= L(i,k) * L(i,k);
		
		L(i,i) = sqrt(tmp);
    }
	
	return(L);
}




//サンプリング関連

//１次元ガウス分布からサンプリング
double SingleGaussSampler(double mu,double sig)
{
	double x1 = double(rand()+1.0)/(RAND_MAX +1.0);
	double x2 = double(rand()+1.0)/(RAND_MAX +1.0);
	
	double y = sqrt(-2*log(x1))*cos(2.0*M_PI*x2) ;
	
	return(mu + sig * y);
}
//1次元逆ガンマ分布からサンプリング
double SingleInvGammaSampler(double kappa,double lambda)
{
	return(1.0 / gengam(kappa,lambda));
}


//平均ベクトルと分散共分散行列に基づく多次元ガウス分布からサンプリング
dcovector MultiGaussSampler(dcovector Mu,dgematrix Sig)
{
	double r1,r2;
	dcovector value(Mu.l);
	dgematrix L = Sig;
	double tmp=0;
	
	for(int i=0;i<Mu.l;i++)
    {
		r1 = double(rand()+1.0)/(RAND_MAX+1.0);
		r2 = double(rand()+1.0)/(RAND_MAX+1.0);
		
		value(i) = sqrt(-2 * log(r1))*cos(2.0*PI*r2);
    }
	
	L.zero();
	for(int i=0;i<Mu.l;i++)
    {
		//i>jの計算
		for(int j=0;j<i;j++)
        {
			tmp = Sig(i,j);
			for(int k=0;k<j;k++)
				tmp -= L(i,k) * L(j,k);
			
			L(i,j) = tmp / L(j,j);
        }
		//i=jの計算
		tmp = Sig(i,i);
		for(int k=0;k<i;k++)
			tmp -= L(i,k) * L(i,k);
		
		L(i,i) = sqrt(tmp);
    }
	
	value = (L*value + Mu);
	
	return(value);
}

//確率を要素とする多項分布からサンプリング
int MultiNominalSampler(dcovector v)
{
	dcovector pr = v;
	double t = 0.0;
	double sum = 0;
	
	for(int i=0;i<v.l;i++)
    {
		sum += v(i);
    }
	pr *= (1.0/sum);
	t = double(rand())/(double(RAND_MAX) + 1.0);
	
	int k =0;
	for(int i=0;i<v.l;i++){
		if(t>=0){
			t=t-pr(i);
			if(t<0){k=i;}
		}
	}
	return k;
}

//Dirichlet分布からサンプリング
dcovector DirichletSampler(dcovector vec)
{
	dcovector theta(vec.l);theta.zero();
	double sum=0;
	
	for(int i=0;i<vec.l;i++)
    {
		if(vec(i) != 0)
			theta(i) = gengam(1.0,vec(i));
		else
			theta(i) = 0;
		sum += theta(i);
    }
	
	for(int i=0;i<vec.l;i++)
		theta(i) = theta(i) / sum;
	
	return(theta);
}

//2項分布からサンプリング
int BinominalSampler(int n,double p)
{
	dcovector x(2);
	x(0) = p;
	x(1) = 1.0 - p;
	
	int a=0;
	for(int i=0;i<n;i++)
		if(0 == MultiNominalSampler(x))
			a++;
	
	return(a);
}

//逆Wishart分布からサンプリング
//パラメータは自由度n,精度行列S
//入力を分散投入に変更2010/10/06taniguchi
//実際には間違っていなかったので，戻す．2011/3/31
dgematrix IWishartSampler(double n,dgematrix S)
{
  S = CPPL::i(S); //precision to covariance
	dgematrix z(S.m,S.m);z.zero();
	dcovector c(S.m);c.zero();
	dgematrix R(S.m,S.m);R.zero();
	
	for(int i=0;i<S.m;i++)
		for(int j=0;j<S.m;j++)
			z(i,j) = gennor(0,1);
	
	for(int i=0;i<S.m;i++)
		c(i) = genchi(n-(double)i);
	
	for(int i=0;i<S.m;i++)
		for(int j=i;j<S.m;j++)
		{
			if(i == j)
				R(j,i) = pow(c(i),0.5);
			else
				R(j,i) = z(j,i);
		}
	dgematrix X = R*CPPL::t(R);
	dgematrix C = cholesky(S);
	dgematrix D = CPPL::t(C) * X * C;
	return (CPPL::i(D));
}

/*
dgematrix IWSample(double n, dgematrix A){
  return (IWishartSampler(n,i(A)));
}
*/

void dge_resize(dgematrix *X, int m, int n,int number ){
	for(int i =0;i<number ;i++){
		X[i].resize(m,n);
		X[i].identity();
	}
}



drovector sum_to_dro(dgematrix X){
	drovector temp(X.n);
	temp.zero();
	for(int i=0;i<X.m;i++){
		temp+=rovec_read(X,i);
	}
	
	return temp;
}


//NBGauss // Gaussian

void NBGauss::resize(int k){
	Mu.resize(k);Mu.zero();
	Sig.resize(k,k);Sig.identity();
	hp_m.resize(k);hp_m.zero();
	hp_S.resize(k,k);hp_S.identity();
	hp_n = 1;
	hp_A.resize(k,k);hp_A.identity();
}

dcovector NBGauss::Sampler(){
	return MultiGaussSampler(Mu,Sig);
	
}

double NBGauss::Probability(dcovector x){
	return Cal_MultiNormLikely(x,Mu,Sig);
}



dcovector NBGauss::UpdateMu(dgematrix X){
	
	dgematrix Sig_bar = i(i(hp_S) + X.m*i(Sig));
	dcovector musum = (i(hp_S)*hp_m + i(Sig)*t(sum_to_dro(X)));
	
	dcovector Mu_bar = Sig_bar*musum;
	
	Mu = MultiGaussSampler(Mu_bar,Sig_bar);
	return Mu;
}

dgematrix NBGauss::UpdateSig(dgematrix X){
	
	double nu = hp_n+X.m;
	
	dgematrix mean_matrix(X.m,X.n);
	
	for (int i=0;i<X.m; i++) {
		drovector tMu = t(Mu);
		vec_set(mean_matrix,i,tMu);
	}
	
	//cout << mean_matrix << endl;
	
	dgematrix Y  = X - mean_matrix; 
	
	dgematrix Delta = hp_A + t(Y)*Y;
	
	Sig = IWishartSampler(nu,Delta);
	return Sig;
}



/// NBMulit // 多項分布

void NBMulti::resize(int k){
	Mu.resize(k);all_set(Mu,1.0/double(k));
	hp_alpha.resize(k);all_set(hp_alpha,DEFAULT_HP_ALPHA);
}

double NBMulti::Probability(int i){
	return Mu(i);
}

int NBMulti::Sampler(){
	return MultiNominalSampler(Mu);
}

dcovector NBMulti::UpdateMu(dcovector count_vec){
	
	Mu = DirichletSampler(count_vec + hp_alpha);
	return Mu;
}


//index とkeyが一致している行だけをとりだして，make matrix 
dgematrix ExtractMatrixByIndex(dgematrix X, vector<int> index, int key){
	
	int count=0;
	
	for(int i=0;i<X.m;i++){
		if(index[i]==key){
			count++;
		}
	}
	
	dgematrix Y(count,X.n);
	
	count=0;
	for(int i=0;i<X.m;i++){
		if(index[i]==key){
			vec_set(Y,count,rovec_read(X,i));
			count++;
		}
	}
	
	
	return Y;
}





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


dgematrix diag(dcovector x){
	dgematrix mat(x.l,x.l);
	
	mat.identity();
	for (int i=0; i<x.l; i++) {
		mat(i,i) = x(i);
	}
	
	return mat;
	
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
	
	
};


dcovector direct_sum(dcovector x , double y){
	dcovector v;v.resize(x.l+1);
	for (int i=0; i<x.l; i++) {
		v(i) = x(i);
	}
	v(x.l) = y;
	
	return v;
	
}

double sum(dcovector x){
	double a = 0;
	for (int i=0; i<x.l; i++) {
		a+= x(i);
	}
	return a;
}

// additional column is introduced
dgematrix ForwardFiltering(NBHmm H, dgematrix X){
	
	int num = H.M.size();
	
	double log_c_ = 0; 
	
	dgematrix F(X.m,num+1); // the last column is for baseline
	F.zero();
	
	dcovector current(num);
	drovector x = rovec_read(X,0);	
	for (int i=0; i<num; i++) {
		current(i) = H.G[i].Probability(t(x));
	}
	
	
	vec_set(F,0,t(direct_sum(current,1)));
	
	H.TM();
	for(int j=1;j<X.m;j++){
		current = H.TM_buffer*current;
		x = rovec_read(X,j);
		for (int i=0; i<current.l; i++) {
			current(i) = H.G[i].Probability(t(x))*current(i);
		}
		double a = sum(current);
		current = (1.0/a)*current;
		log_c_ = log_c_ + log10(a);//similar to hamahata
		vec_set(F,j,t(direct_sum(current,log_c_)));
	}
	
	return F;
}



dcovector subvector(dcovector x, int l){
	dcovector v(l);
	for (int i=0; i<l; i++) {
		v(i) = x(i);
	}
	return v;
}
drovector subvector(drovector x, int l){
	drovector v(l);
	for (int i=0; i<l; i++) {
		v(i) = x(i);
	}
	return v;
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


dcovector dco(vector<int> v){
	dcovector x(v.size());
	for (int i=0; i<x.l; i++) {
		x(i) = v[i];
	}
	return x;
	
}

dgematrix TransitionCount(vector<int> s, int states){
	dgematrix mat(states,states);
	mat.zero();
	for (int i=0; i< s.size()-1; i++) {
		mat(s[i+1],s[i])++;
	}
	return mat;
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





