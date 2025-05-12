#include <numeric>
#include <stdio.h>
#include <cblas.h>
#include <iostream>
using namespace std;

typedef vector<int64_t> polynomial;

// polynomial a mod q
//a[i] in [-q/2,q/2)
inline void modq_poly(polynomial &a, int n, int64_t q) {
  
  // INPUT: polynomial a, modulus q
  // OUTPUT: (modified) a mod q, stored in a

  for (int i = 0; i < n; ++i) {

    while (a[i]<0) {
      a[i] += q;
    }

    while (a[i]>=q) {
      a[i] -= q;
    }

    if(a[i] > (q-1)/2){
      a[i] -=q;
    }

  }
  
}

inline void modq_poly_large(polynomial &a, int n, int64_t q) {
  
  // INPUT: polynomial a, modulus q
  // OUTPUT: (modified) a mod q, stored in a

  for (int i = 0; i < n; ++i) {

    if(a[i]<0){
      int64_t temp = -1 * a[i];
      a[i] += q*(temp/q+1);
    }

    if(a[i] >= q){
      a[i] -= q*(a[i]/q);
    }

    if(a[i] > (q-1)/2){
      a[i] -=q;
    }

  }
  
}

// Addition between polynomials
 extern void add_poly(polynomial &a, const polynomial & b, int n, int64_t q) {
  // INPUT: polynomials a and b, modulus q
  // OUTPUT: (a + b) mod q, stored in a
  for(int i = 0 ; i < n ; ++i){
    a[i] += b[i];
  }
  // mod q
  modq_poly_large(a,n, q);

}

 // Scaling
 extern void multi_scale_poly(int64_t t, polynomial &a, int n, int64_t q) {
  // INPUT: scaler t, polynomial a, modulus q
  // OUTPUT: t*a mod q, stored in a

  for (int i = 0; i < n; ++i) {
    a[i] *= t;
  }
  modq_poly_large(a,n, q);
}

extern polynomial mul_poly( const polynomial& aa, const polynomial& bb,int n, int64_t q) {
  polynomial c(n,0);
  for(int i = 0 ; i <n ; ++i){
    for(int j = 0 ; j < n ; ++j){
      if(i+j < n){
        c[i+j] += aa[i]*bb[j];
      }
      else{
        c[i+j-n] -= aa[i]*bb[j];
      }
    }
  }
  modq_poly_large(c,n,q);
  return c;
}

extern vector<double> mul_poly_double( const polynomial& aa, const polynomial& bb,int n) {
  vector<double> c(n,0);
  for(int i = 0 ; i <n ; ++i){
    for(int j = 0 ; j < n ; ++j){
      if(i+j < n){
        c[i+j] += (double)aa[i]*(double)bb[j];
      }
      else{
        c[i+j-n] -= (double)aa[i]*(double)bb[j];
      }
    }
  }
  //modq_poly(c,n,q);
  return c;
}

// RLWE key generation algorithm
// INPUT: dimension n
// OUTPUT: a degree-(n-1) polynomial s, with coeffients from Uniform Random Distribution on {-1,0,1}
extern polynomial RLWE64_KeyGen(int n) {
  unsigned int len_out = n;
  unsigned char seed[SEED_LEN];
  int *x = new int [n];
  random_bytes(seed,SEED_LEN);
  int r = gen_ternary(x,len_out,seed);
  if (r == 1){
    cout <<"Error in RLWE Key Generation: NULL "<<endl;
  }
  if(r == 2){
    cout <<"Error in RLWE Key Generation: Length" <<endl;
  }

  polynomial s(n,0);
  //s[0] = (int64_t)n - 1;
  for (int i = 0; i < n; ++i) {
    s[i] = x[i];
  }
  return s;
}

// RLWE encryption algorithm
// INPUT: dimension n, modulus q, variance 2^(-k), k >= 1, plaintext (polynomial) m, RLWE key s
// OUTPUT: RLWE ciphertext ct = (a, b) = (a, alpha*m + e - a * s)
extern vector<polynomial> RLWE64_Enc(int n, int64_t q, double alpha, double k, const vector<double> & m, const polynomial & s) {

  vector<polynomial> ct;

  //generate random a
  double logq=log(q)/log(2.0);
  int int_logq = (int) logq;
  if(logq > (double)int_logq){
    int_logq++;
  }

  unsigned int len_out = n;
  unsigned char seed[SEED_LEN];
  int64_t *array_a = new int64_t [len_out];
  random_bytes(seed,SEED_LEN);
  int r = gen_uniform_int64(array_a, len_out, q, int_logq, seed);

  if (r == 1){
    cout <<"Error in generation random array a of RLWE Encryption: NULL "<<endl;
  }
  if(r == 2){
    cout <<"Error in generation random array a of RLWE encryption: modulus" <<endl;
  }

  polynomial a(n,0);
  //a[0] = (int64_t)n - 1;
  for (int i = 0; i < n; ++i) {
     a[i]=array_a[i]-q/2;
  }
  
  ct.push_back(a);

  //compute - a * s
  polynomial as = mul_poly(s,a,n,q);
  //print_polynomial(as);
  //cout <<"as computed."<<endl;
  multi_scale_poly(-1, as,n,q);

  //cout <<"as computed."<<endl;

  //generate  error e
  //generate error array e
  int *array_e = new int [len_out];
  random_bytes(seed,SEED_LEN);
  r = gen_discrete_normal(array_e, len_out, k, seed);
  if (r == 1){
    cout <<"Error in generation a random error of LWE Encryption: NULL "<<endl;
  }
  if(r == 2){
    cout <<"Error in generation a random error of LWE Encryption: modulus" <<endl;
  }

  polynomial as_m(n);
  for (int i = 0; i < n; ++i) {
    //as_m[i] = as[i] + alpha*m[i];
    if(array_e[i]>10 || array_e[i]<-10){
      cout <<"error > 10"<<endl;
    }
    double tempm = alpha*m[i];
    int64_t r_tempm = (int64_t)tempm;
    if(tempm > 0 && tempm - (double)r_tempm >= 0.5){
      r_tempm ++;
    }
    else if(tempm < 0 && tempm - (double)r_tempm <= -0.5){
      r_tempm --;
    }
    as_m[i] = as[i] + r_tempm;
    //as_m[i] = as[i] + r_tempm+array_e[i];
  }
  modq_poly_large(as_m, n,q);
 // cout <<"as+m+e computed."<<endl;
  ct.push_back(as_m);

  return ct;
}

// RLWE Decryption algorithm  
// INPUT: dimension n, modulus q, RLWE key s, RLWE ciphertext ct (a, b)
// OUTPUT: polynomial (1/alpha)(b + a * s)
extern vector<double> RLWE64_Dec(int n, int64_t q, double alpha, const polynomial & s, const vector<polynomial> & ct) {
  //compute as
  polynomial as = mul_poly(s, ct[0],n,q);
  //compute b+as
  add_poly(as, ct[1],n,q);

  vector<double> output(n);

  //(1/alpha)*as
  for(int i = 0 ; i < n ; ++i){
    output[i] = (1.0/alpha)*(double)as[i] ;
  }

  return output;
}

// modulus Switching, mod q1 -> mod q2  
// INPUT: first modulus q1, second modulus q2, dimension n, LWE ciphertext ct with modulus q1
// OUTPUT: LWE ciphertext with modulus q2
vector<polynomial> RLWE64_Rounding(int64_t q1, int64_t q2, int n, const vector<polynomial>& ct) {
  vector<polynomial> ans=ct;
  for (int i = 0; i < 2; ++i) {
    for (int j = 0 ; j < n ; ++j){
      if(ans[i][j] < 0){
        ans[i][j] += q1;
      }
      __int128 temp = (__int128)ans[i][j]*(__int128)q2 + (__int128)q1/2;
      temp /= (__int128)q1;
      ans[i][j] = (int64_t)temp;
    }
  }
  return ans;
}

//U=t*d, M = d*n, output = t*n = t ciphertexts
extern vector<vector<polynomial>> RLWE64_matrix_mul(int t, int d, int n, 
  int64_t q, const vector<vector<polynomial>>& ct, 
  const vector<vector<float>> &matrixU){
  //construct A,B,U
  //float A[d*n];
  float *A = new float[d*n];
  for (int i = 0; i < d; ++i){
    for (int j = 0 ; j < n ; ++j){
      A[i*n+j] = (float)ct[i][0][j];
    }
  }
 // cout <<"A"<<endl;

  

  //float U[t*d];
  float *U = new float[t*d];
  for (int i = 0; i < t; ++i){
    for (int j = 0 ; j < d ; ++j){
      U[i*d+j] = (float)matrixU[i][j];
    }
  }

 // cout <<"U"<<endl;

  //compute UA, UB
  //float UA[t*n];
  float *UA = new float[t*n];
  cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, t, n, d, 1, 
    U, d, A, n, 0, UA, n); 

  delete []A;

  //generate t ct from UA and UB
  vector<vector<polynomial>> output(t, vector<polynomial>(2,polynomial(n,0)));
  for (int i = 0; i < t; ++i){
    for (int j = 0 ; j < n ; ++j){
      output[i][0][j] = (int64_t)UA[i*n+j];
    }
    modq_poly_large(output[i][0],n,q);
  }

  delete []UA;


  //float B[d*n];
  float *B = new float[d*n];
  for (int i = 0; i < d; ++i){
    for (int j = 0 ; j < n ; ++j){
      B[i*n+j] = (float)ct[i][1][j];
    }
  }
//  cout <<"B"<<endl;

  //float UB[t*n];
  float *UB = new float[t*n];
  cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, t, n, d, 1, 
    U, d, B, n, 0, UB, n);

  delete []B;
  
  for (int i = 0; i < t; ++i){
    for (int j = 0 ; j < n ; ++j){
      output[i][1][j] = (int64_t)UB[i*n+j];
    }
    modq_poly_large(output[i][1],n,q);
  }

  return output;



}