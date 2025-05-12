using namespace std;
using namespace seal;

inline
long power(long x, long y, long m)
{
    if (y == 0)
        return 1;
    long p = power(x, y / 2, m) % m;
    p = (p * p) % m;
 
    return (y % 2 == 0) ? p : (x * p) % m;
}

inline
long modInverse(long a, long m)
{
    return power(a, m - 2, m);
}

vector<vector<int>> generateMatrixU_transpose(int n, const int q = 65537) {
    vector<vector<int>> U(n,  vector<int>(n));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == 0) {
                U[i][j] = (int) power(3, j, q);
            } else if (i == n/2) {
                U[i][j] = (int) modInverse(U[0][j], q);
            } else {
                U[i][j] = (int) power(U[i-1][j], 3, q);
            }
        }
    }
    return U;
}

vector<vector<uint64_t>> cal_dcd_matrix(int degree, const SEALContext& seal_context){

  BatchEncoder batch_encoder(seal_context);

  vector<vector<uint64_t>> dcdmatrix(degree, vector<uint64_t>(degree,0));

  for(int i = 0 ; i < degree ; ++i){
    Plaintext plainInd;
    plainInd.resize(degree);
    plainInd.parms_id() = parms_id_zero;
    for (int j = 0; j < (int) degree; ++j) {
        plainInd.data()[j] = 0;
    }
    plainInd.data()[i] = 1;

    vector<uint64_t> podmatrix;
    batch_encoder.decode(plainInd,podmatrix);

    for(int j = 0 ; j < degree ; ++j){
      dcdmatrix[j][i] = podmatrix[j];
    }
    

  }

  return dcdmatrix;

  
}


vector<Plaintext> ecd_dcd_matrix(int degree, const vector<vector<uint64_t>>& dcdmatrix, 
  const SEALContext& seal_context, const Ciphertext & RLWE_ct){

  Evaluator evaluator(seal_context);
  BatchEncoder batch_encoder(seal_context);
  size_t slot_count = batch_encoder.slot_count();

  int sq_degree_div_2 = sqrt(degree/2);
  
  vector<Plaintext> output(degree);

  #pragma omp parallel for 

  for (int i = 0; i < sq_degree_div_2; ++i){
    for (int j = 0 ; j < 2*sq_degree_div_2 ; ++j){

      vector<uint64_t> temp(degree);

      for (int k = 0; k < degree; ++k){
        int index_row = (k - i)%(degree/2);
        if(index_row < 0){
          index_row += degree/2;
        }
        if(k >= degree/2){
          index_row += degree/2;
        }

        int index_col = (k+j*sq_degree_div_2)%(degree/2);
        if(j < sq_degree_div_2){
          if(k >= degree/2){
            index_col += degree/2;
          }
        }
        else{
          if(k < degree/2){
            index_col += degree/2;
          }
        }
        temp[k] = (uint64_t)dcdmatrix[index_row][index_col];
      }

     // if(j == 0) cout <<temp[10]<<" ";

      Plaintext temp_pt;
      batch_encoder.encode(temp,temp_pt);
      evaluator.transform_to_ntt_inplace(temp_pt, RLWE_ct.parms_id());
      output[i * 2*sq_degree_div_2 + j] = temp_pt;

    }
    
  }

  return output;

}

extern Ciphertext Slot_to_coeff(const Ciphertext & RLWE_ct, const vector<Plaintext> & M, 
const GaloisKeys & RotK, const SEALContext& seal_context, const SecretKey & sk){

  Evaluator evaluator(seal_context);
  int degree = M.size();
  int sq_degree_div_2 = sqrt(degree/2);

  struct timeval tstart1, tend1;

  gettimeofday(&tstart1,NULL);

  //baby step 1

  vector<Ciphertext> c_g(sq_degree_div_2,RLWE_ct);
  evaluator.transform_to_ntt_inplace(c_g[0]);   

  #pragma omp parallel for         

  for(int i = 1 ; i < sq_degree_div_2 ; ++i){
    evaluator.rotate_rows_inplace(c_g[i], sq_degree_div_2 * i, RotK);
    evaluator.transform_to_ntt_inplace(c_g[i]);
  }

  //baby step 2

  Ciphertext RLWE_ct_inv(RLWE_ct);
  evaluator.rotate_columns_inplace(RLWE_ct_inv,RotK);

  vector<Ciphertext> c_g_2(sq_degree_div_2,RLWE_ct_inv);
  evaluator.transform_to_ntt_inplace(c_g_2[0]);   

  #pragma omp parallel for           

  for(int i = 1 ; i < sq_degree_div_2 ; ++i){
    evaluator.rotate_rows_inplace(c_g_2[i], sq_degree_div_2 * i, RotK);
    evaluator.transform_to_ntt_inplace(c_g_2[i]);
  }

  gettimeofday(&tend1,NULL);
  double baby_time = tend1.tv_sec-tstart1.tv_sec+(tend1.tv_usec-tstart1.tv_usec)/1000000.0;
  //cout <<"Baby step time = "<<baby_time<<endl;

  //gaint step

  gettimeofday(&tstart1,NULL);

  vector<Ciphertext> output(sq_degree_div_2);

  #pragma omp parallel for 
  
  for (int i = 0; i < sq_degree_div_2; ++i){
    for (int j = 0; j < sq_degree_div_2; ++j){
      if(j == 0){
        evaluator.multiply_plain(c_g[j],M[i*2*sq_degree_div_2+j],output[i]);
      }
      else{
        Ciphertext temp;
        evaluator.multiply_plain(c_g[j],M[i*2*sq_degree_div_2+j],temp);
        evaluator.add_inplace(output[i], temp);
      }
    }
    for (int j = 0; j < sq_degree_div_2; ++j){
      Ciphertext temp;
      evaluator.multiply_plain(c_g_2[j],M[i*2*sq_degree_div_2+j+sq_degree_div_2],temp);
      evaluator.add_inplace(output[i], temp);
    }
  }

  for (int i = 0; i < sq_degree_div_2; i++) {
    evaluator.transform_from_ntt_inplace(output[i]);
  }

  for (int i = sq_degree_div_2-1; i > 0; --i) {
    evaluator.rotate_rows_inplace(output[i], 1, RotK);
    evaluator.add_inplace(output[i-1], output[i]);
  }

  gettimeofday(&tend1,NULL);
  double gaint_time = tend1.tv_sec-tstart1.tv_sec+(tend1.tv_usec-tstart1.tv_usec)/1000000.0;
  //cout <<"Gaint step time = "<<gaint_time<<endl;

  return output[0];

}