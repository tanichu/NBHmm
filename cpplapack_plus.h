/*#pragma once
#pragma comment(lib,"libF77.lib")
#pragma comment(lib,"libI77.lib")
#pragma comment(lib,"blas.lib")
#pragma comment(lib,"clapack.lib")
*/
	
#define CSV_STRING vector< vector< string > >
#define _USE_MATH_DEFINES 
	
#define PI 3.141592	
#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>

#include <cpplapack.h>//cpplapack
#include <randlib.h> // randlib


#include <time.h>//時間関数利用　乱数初期化など
//#include <conio.h>//コンソールで一時止めるgetch()を利用するため
#include <math.h>//数学関数
#include <string.h>//文字列処理用のクラスstring利用のため
#include <fstream>//ファイルストリーム
using namespace std;
using namespace CPPL;

//自作の文字列変換用
std::string itos(int n); // int -> string
std::string ftos(double n);// float,double -> string
string TrimRight( string arg ); //string の右側のスペースなどをとる．
string TrimLeft( string arg );//string の左側のスペースなどをとる．
string Trim( string arg);//string の両側のスペースなどをとる．

string ScanTillRet(string str); // Return迄の文字列を取得
string ScanTillSpace(string str); // space迄の文字列を取得
string ScanTillTab(string str); // tab迄の文字列を取得
string ScanTillComma(string str);// comma迄の文字列を取得
string DelTillRet(string str);// Return迄の文字列を消去
string DelTillSpace(string str);// space迄の文字列を消去
string DelTillTab(string str);// Return迄の文字列を消去
string DelTillComma(string str);// Return迄の文字列を消去

//他にもデフォルトで
// atof // string -> float
// atoi // string -> int
//が存在します．数字と文字列の相互変換に使ってください．

//void DEBUG_TRACE(string x);
vector<string> split(const string& str, const string& delimiter);

vector< vector<string> > csv_reader(const char *filename);
vector< vector<string> > ssv_reader(const char *filename);




double gauss_rand(); //分散１平均０のガウス乱数を発生
double gauss_rand(double ave,double sig);//平均ave,標準偏差sigのガウス乱数を発生
dcovector gauss_rand(dcovector ave,double sig);//ベクトルで平均ave,全次元独立で偏差sigのガウス乱数
dcovector gauss_rand(dcovector ave,dcovector sig);//ベクトルで平均ave,全次元独立で各次元偏差sigのガウス乱数
dgematrix gauss_rand(dgematrix ave,double sig);//行列の各成分に偏差sigの乱数を載せる

double LGF(dcovector x,dcovector m,dcovector s);//対数化したガウス関数，中心m，各次元独立で標準偏差sのガウス関数をlogに入れたもの．ガウス関数のかけ算をやる際などに，対数空間内で和算で代用したほうが精度が出る．
double LGF(dcovector x, dcovector m,dgematrix q); //
double LGF(dcovector x, dcovector m,dgematrix iq,double detq);


dcovector ng_calculation(dcovector *m,dcovector *s,int N,dcovector x);//複数のガウス分布の中心，標準偏差，ガウス分布の数Nを渡すことで正規化ガウス関数(NGnet)の出力を出す．出力は各基底からの出力値
dcovector rbf_calculation(dcovector *m,dcovector *s,int N,dcovector x);//上のＲＢＦ版

void restrict(dcovector &v, double min,double max);//ベクトルの値をmin,maxの間に制限する（外れ値処理などに）
void restrict(drovector &v, double min,double max);//ベクトルの値をmin,maxの間に制限する
void restrict(dgematrix& A,double min,double max);//行列の値をmin,maxの間に制限する
void  each_square(dcovector &v);//各次元について自乗する．
void  each_abs(dcovector &v);//各次元について絶対値を与える．
void scatter();//乱数初期化　時間関数を呼び出してrand関数を初期化する．プログラムの一番始めでとりあえず呼ぶべし．

double frobnorm(dgematrix A);//行列のフロベニウスノルムを計算
double trace(dgematrix A);//行列のトレースを計算

void read_multi_vector(dcovector *x,int num, int dim, char *filename);//ファイルから複数のベクトルを読み取る．
void read_multi_vector(drovector *x,int num, int dim, char *filename);//ファイルから複数のベクトルを読み取る．
void read_multi_matrix(char *filename,dgematrix *A,int n);//ファイルから複数の行列を読み取る．


dgematrix product(dcovector c,drovector r);//列ベクトルと行ベクトルの積で行列を造る．
//dgematrix sup::operator*(dcovector c,drovector r);

void vec_set(dgematrix &A, int k , dcovector x); // 行列の一部にベクトルをセット

void vec_set(dgematrix &A, int k , drovector x);
dcovector covec_read(dgematrix A, int k);
drovector rovec_read(dgematrix A, int k);

int ProbSelect(dcovector v);//各次元の値を確率と見なし，その確率に沿ってサイコロを振り次元を選択する．

double max_value(dcovector v);//各次元値の中のの最大値を返す．
double min_value(dcovector v);//各次元値の中の最小値を返す
int argmax(dcovector v);//最大値を持つ次元を返す
int argmin(dcovector v);//最小値を持つ次元を返す
dcovector argmax_vector(dcovector v);//
dcovector argmin_vector(dcovector v);
dcovector e_greedy(dcovector a,double greedy_ratio);//ProbSlectのe-greedy版，最大値を持つ次元を1-(N-1)eの確率，他をe=greedy_ratioの確率で選択．


dcovector normalize(dcovector x);//xを各次元の和が１になるように正規化する．
dcovector each_product(dcovector a,dcovector b);//
double det(dgematrix q); //行列式を返す
double gaussian(dcovector x, dcovector m,dgematrix q); //
double gaussian(dcovector x, dcovector m,dgematrix iq,double detq);

double norm2(dcovector x);//二乗ノルムを返す


//dgematrix inv(dgematrix A);

string csv(dcovector x);//csvの形式の文字列に出力する． fout << csv(X) <<endl;といった使い方
string csv(drovector x);//同上
string csv(dgematrix x);//同上

void self_add(dgematrix &A,dgematrix B);//selfシリーズは演算子を用いず，自らに計算結果を代入する．
void self_add(dcovector &x,dcovector y);//多少の速度向上は期待できるかもしれないが
void self_set(dgematrix &A,dgematrix B);//そんなに利用価値はない．
void self_set(dcovector &x,dcovector y);//


void all_set(dgematrix &x,double value);//行列の全成分をvalueにセット
void all_set(dcovector &x,double value);//行列の全成分のvalueにセット

dgematrix submatrix(dgematrix A, int L, int R, int T, int B);

dgematrix lean_on (dgematrix A, dgematrix B);
dgematrix put_on (dgematrix A, dgematrix B);

//dgematrix operator = (dcovector v);



///////routines used in Hamahata's sticky-HDP-HMM program

//--------------------------------------------------------
//-------------------------関数----------------------------
//--------------------------------------------------------
//行列にベクトルを挿入
dgematrix Insert_mat(dgematrix mat,int col,dcovector vec);
//行列からベクトルを抽出
dcovector Extract_vec(dgematrix mat,int col);
//ベクトルのノルムを計算
double Cal_norm(dcovector vec);
//行列にベクトルを挿入して足し合わせる
dgematrix InsertAdd_VecToMat(dgematrix mat,int col,dcovector vec);

//正則行列かどうかのチェック（正則なら0を返す）
int check_Regularization(dgematrix mat);


//多次元正規分布の尤度を計算
double Cal_MultiNormLikely(dcovector x,dcovector mu,dgematrix sig);
//コレスキー分解
dgematrix cholesky(dgematrix mat);



//サンプリング関連
//多次元ガウス分布からサンプリング
dcovector MultiGaussSampler(dcovector Mu,dgematrix Sig);
//1次元ガウス分布からサンプリング
double SingleGaussSampler(double mu,double sig);
//1次元逆ガンマ分布からサンプリング
double SingleInvGammaSampler(double kappa,double lambda);
//多項分布からサンプリング
int MultiNominalSampler(dcovector v);
//ディリクレ分布からサンプリング
dcovector DirichletSampler(dcovector vec);
//2項分布からサンプリング
int BinominalSampler(int n,double p);
//逆Wishart分布からサンプリング
dgematrix IWishartSampler(double n,dgematrix S);



///NBHmm 開発時に導入
dcovector dco(vector<int> v);
void dge_resize(dgematrix *X, int m, int n,int number );

drovector sum_to_dro(dgematrix X);

dgematrix diag(dcovector x);
	
dgematrix ExtractMatrixByIndex(dgematrix X, vector<int> index, int key);

dcovector direct_sum(dcovector x , double y);
double sum(dcovector x);
dgematrix TransitionCount(vector<int> s, int states);


dcovector subvector(dcovector x, int l);
drovector subvector(drovector x, int l);




	

#define DEFAULT_HP_ALPHA 0.1 
class NBMulti{
	
public:
	dcovector Mu;
	
	dcovector UpdateMu(dcovector count);
	int Sampler();
	void resize(int k);
	
	double Probability(int i);
	//Hyperparameters                                                             
	dcovector hp_alpha;
	
};

	

	
class NBGauss{
	
public:
	dcovector Mu;
	dgematrix Sig;
	
	dcovector UpdateMu(dgematrix X);
	dgematrix UpdateSig(dgematrix X);
	dcovector Sampler();
	
	double Probability(dcovector x); 
	
	void resize(int k);
	//Hyperparameters                                                             
	//for mean                                                                   
	dcovector hp_m;
	dgematrix hp_S;
	// for variance                                                               
	double hp_n;
	dgematrix hp_A;
	
};

//NBHmm  Forwardfiltering - Backward sampling のBayseHMM
class NBHmm  {
public:
	vector<NBGauss> G;//Gaussian distribution
	vector<NBMulti> M;// Transition multinomial distribution
	
	void resize(int num_states, int dim_output);
	dgematrix TM_buffer; // used to buffering TM for reducing repetedly estimation of TM.
	dgematrix TM();// shows transition matrix
	
	
	void Update(dgematrix Y,vector<int> label);
	
	void read_Mu(const char *filename);
	void read_diag_Sig(const char *filename);
	
	void read_TM(const char *filename);// read TM from file
};





dgematrix ForwardFiltering(NBHmm H, dgematrix X);
vector<int> BackwardSampling(NBHmm H,dgematrix F);
vector<int> GenerateStates(NBHmm H, int length, int initial_state);
dgematrix GenerateObservations(NBHmm H, vector<int> s);


