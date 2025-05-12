using namespace std;
using namespace seal;

inline vector<uint64_t> Rot_Vec(vector<uint64_t> input, int k){
  int len = input.size();
  vector<uint64_t> out(len);
  for(int j = 0 ; j <2 ; ++j){
    for(int i = 0 ; i < len/2 ; ++i){
      int index = i+k;
      if(index >= len/2){
        index -= len/2;
      }
      out[index+j*len/2] = input[i+j*len/2];
    }
  }
  return out;
}


extern vector<Plaintext> ecd_M(const vector<vector<uint64_t>> & M, int N, const SEALContext& seal_context){
  
  BatchEncoder batch_encoder(seal_context);
  size_t slot_count = batch_encoder.slot_count();

  int M_l = M.size(); 
  int M_n = M[0].size();

  int n1,n2;
  if(M_l >= M_n){
    n1 = M_l;
    n2 = M_n;
  }
  else{
    n1 = M_n;
    n2 = M_l;
  }

  //construct n2 vectors
  vector<vector<uint64_t>> M2(n2, vector<uint64_t>(n1,1));

  #pragma omp parallel for               

  for(int i = 0 ; i < n2 ; ++i){
    for(int j = 0 ; j < n1 ; ++j){
      //index1 = j mod M_l
      int index1 = j;
      while(index1 >= M_l){
        index1 -= M_l;
      }
      //index2 = j+i mod M_n
      int index2 = j+i;
      while(index2 >= M_n){
        index2 -= M_n;
      }

      M2[i][j] = M[index1][index2];
    }
  }

   double sqrt_n2 = sqrt((double)n2);

  //rounding sqrt_n2
  int r_sqrt_n2 = (int)sqrt_n2;
  if(sqrt_n2-(double)r_sqrt_n2 > 0){
    r_sqrt_n2 ++;
  }

  //giant step
  double b = (double)n2/(double)r_sqrt_n2;
  int r_b = (int)b;
  if(b-(double)r_b > 0){
    r_b++;
  }

  vector<Plaintext> out(N);

  for(int i = 0 ; i < r_sqrt_n2;++i){
    batch_encoder.encode(M2[i], out[i]);

  }

  #pragma omp parallel for                 

   for(int j = 1; j < r_b-1 ; ++j){
    for (int i = 0; i < r_sqrt_n2; ++i){
      int index = j*r_sqrt_n2+i;
      int index2 = j*r_sqrt_n2;
      while(index2 >= M_n/2){
        index2 -= M_n/2;
      }
      batch_encoder.encode(Rot_Vec(M2[index],index2), out[index]);
    }
  }

  for (int i = 0; i < r_sqrt_n2; ++i){
      int index = (r_b-1)*r_sqrt_n2+i;
      if(index >= M_l){
        break;
      }
      int index2 = (r_b-1)*r_sqrt_n2;
      while(index2 >= M_n/2){
        index2 -= M_n/2;
      }
      batch_encoder.encode(Rot_Vec(M2[index],index2), out[index]);
  }


  //check the correctness
  bool cor = true;
  for (int i = 0; i < r_sqrt_n2; ++i){
    vector<uint64_t> dcd_result;
    batch_encoder.decode(out[i],dcd_result);
    for (int j = 0; j < slot_count; ++j){
      if (dcd_result[j] != M2[i][j]){
        cor = false;
        cout <<dcd_result[j] <<" "<<M2[i][j]<<endl;
      }
    }
  }
  cout <<"Checking for M2[0] - M2["<<r_sqrt_n2-1<<"]........"<<cor<<endl;

  
  return out;

}

extern Ciphertext LT_ecd_M(const Ciphertext & RLWE_ct, const vector<Plaintext> & M, int N,
const GaloisKeys & RotK, const SEALContext& seal_context, const SecretKey & sk){

  Evaluator evaluator(seal_context);
  

  //for test 
  Decryptor decryptor(seal_context, sk);
  BatchEncoder batch_encoder(seal_context);

  int n2 = M.size();
  size_t slot_count = batch_encoder.slot_count();
 // cout <<"size of M = "<<n2<<", slot count = "<<slot_count<<endl;

  double sqrt_n2 = sqrt((double)n2);
 // cout <<"sqrt n2 = "<<sqrt_n2<<endl;

  //rounding sqrt_n2
  int r_sqrt_n2 = (int)sqrt_n2;
  if(sqrt_n2-(double)r_sqrt_n2 > 0){
    r_sqrt_n2 ++;
  }


  //baby step

  vector<Ciphertext> c_g(r_sqrt_n2, RLWE_ct);

  #pragma omp parallel for               

  for(int i = 1 ; i < r_sqrt_n2 ; ++i){
    evaluator.rotate_rows_inplace(c_g[i], i, RotK);
   // evaluator.transform_to_ntt_inplace(c_g[i]);
  }

 // cout <<"Baby step end. "<<endl;

/*
  //test correctness of rotation
  cout <<"Decryption then decoding the rotated ct: "<<endl;
  for (int i = 0; i < r_sqrt_n2; ++i){
    cout <<"i = "<<i<<": ";
    Plaintext plain_result;
    decryptor.decrypt(c_g[i],plain_result);
    vector<uint64_t> pod_matrix;
    batch_encoder.decode(plain_result,pod_matrix);
    for(int j = 0 ; j < 10 ; ++j){
      cout <<pod_matrix[j]<<" ";
    }
    cout <<"...... ";
    for (int j = slot_count-10; j < slot_count; ++j){
      cout <<pod_matrix[j]<<" ";
    }
    cout <<endl;
  }
  */

  //gaint step
  double bb = (double)n2/(double)r_sqrt_n2;
  int r_b = (int)bb;
  if(bb-(double)r_b > 0){
    r_b++;
  }


  Ciphertext out(c_g[0]);
 // evaluator.multiply_plain_inplace(out,M[0]);

  //cout <<0<<endl;
  int non_zero = -1;

  vector<Ciphertext> tempct(r_sqrt_n2,c_g[0]);
  vector<bool> iszero_vec(r_sqrt_n2,true);

  for (int i = 0; i < r_sqrt_n2; ++i){
    //cout <<i<<" ";
    //check whether M[i] is zero vector
   // bool iszero = true;
    for(int k = 0 ; k < slot_count ; ++k){
      if(M[i][k] != 0){
        iszero_vec[i] = false;
        break;
      }
    }
    if(iszero_vec[i] == true){
      //cout <<"M_"<<i<<" is zero vector"<<endl;
      //continue;
    }
    else{
      tempct[i]=c_g[i];
      evaluator.multiply_plain_inplace(tempct[i],M[i]);
      
      
    }
    
  }

  for(int i = 0 ; i < r_sqrt_n2 ; ++i){
    if(iszero_vec[i] == true){
      //cout <<"M_"<<i<<" is zero vector"<<endl;
      continue;
    }
    else{
      if(non_zero == -1){
           out = tempct[i];
           non_zero = 1;
      }
      else{
        evaluator.add_inplace(out,tempct[i]);
      }
    //  evaluator.transform_from_ntt_inplace(out);
    }
  }

  vector<Ciphertext> tempout(r_b-1,c_g[0]);
  vector<int> zero_number(r_b-1,0);

  #pragma omp parallel for

  for(int j = 1; j <r_b ; ++j){
   // cout <<j<<endl;
   // int zero_number = 0;
    int first_non_zero = -1;

    //Ciphertext tempout(c_g[0]);
    //evaluator.multiply_plain_inplace(tempout,M[j*r_sqrt_n2]);

    for(int i = 0; i < r_sqrt_n2 ; ++i){
     // cout <<i<<" ";
      int index = j*r_sqrt_n2+i;
      if(index >= n2){
        break;
      }

      bool iszero = true;
      for(int k = 0 ; k < slot_count ; ++k){
        if(M[j*r_sqrt_n2+i][k] != 0){
          iszero = false;
          break;
        }
      }
      if(iszero == true){
        //cout <<"M_"<<j*r_sqrt_n2+i<<" is zero vector"<<endl;
        zero_number[j-1] ++;
        continue;
      }
      else{
        Ciphertext tempct = c_g[i];
        evaluator.multiply_plain_inplace(tempct,M[j*r_sqrt_n2+i]);
        if(first_non_zero == -1){
           tempout[j-1] = tempct;
           first_non_zero = 1;
        }
        else{
          evaluator.add_inplace(tempout[j-1],tempct);
        }
      }
      
    }

  //  evaluator.transform_from_ntt_inplace(tempout[j-1]);
   // cout <<endl;
    if(zero_number[j-1] == r_sqrt_n2){
      continue;
    }
    else{
      int rot_index = j*r_sqrt_n2;
      while(rot_index >= n2/2){
        rot_index -= n2/2;
      }
      evaluator.rotate_rows_inplace(tempout[j-1],rot_index,RotK);
      
    }
    
  }

  for(int i = 1; i < r_b ; ++i){
    if(zero_number[i-1] == r_sqrt_n2){
      continue;
    }
    else{
      evaluator.add_inplace(out,tempout[i-1]);
    }
  }

  return out;

}


//================== pre stored ntt LT=========
extern vector<Plaintext> ecd_M_ntt(const vector<vector<uint64_t>> & M, int N, 
  const SEALContext& seal_context, const Ciphertext & RLWE_ct){
  
  Evaluator evaluator(seal_context);
  BatchEncoder batch_encoder(seal_context);
  size_t slot_count = batch_encoder.slot_count();

  int M_l = M.size(); 
  int M_n = M[0].size();

  int n1,n2;
  if(M_l >= M_n){
    n1 = M_l;
    n2 = M_n;
  }
  else{
    n1 = M_n;
    n2 = M_l;
  }

  //construct n2 vectors
  vector<vector<uint64_t>> M2(n2, vector<uint64_t>(n1,1));

  //#pragma omp parallel for               

  for(int i = 0 ; i < n2 ; ++i){
    for(int j = 0 ; j < n1 ; ++j){
      //index1 = j mod M_l
      int index1 = j;
      while(index1 >= M_l){
        index1 -= M_l;
      }
      //index2 = j+i mod M_n
      int index2 = j+i;
      while(index2 >= M_n){
        index2 -= M_n;
      }

      M2[i][j] = M[index1][index2];
    }
  }

   double sqrt_n2 = sqrt((double)n2);
 // cout <<"sqrt n2 = "<<sqrt_n2<<endl;

  //rounding sqrt_n2
  int r_sqrt_n2 = (int)sqrt_n2;
  if(sqrt_n2-(double)r_sqrt_n2 > 0){
    r_sqrt_n2 ++;
  }

  //giant step
  double b = (double)n2/(double)r_sqrt_n2;
  int r_b = (int)b;
  if(b-(double)r_b > 0){
    r_b++;
  }

  //cout <<r_sqrt_n2<<" "<<r_b<<endl;

  vector<Plaintext> out(N);

  for(int i = 0 ; i < r_sqrt_n2;++i){
    batch_encoder.encode(M2[i], out[i]);
    evaluator.transform_to_ntt_inplace(out[i], RLWE_ct.parms_id());

  }
//  cout <<"complete 0 - "<<r_sqrt_n2-1<<endl;

  #pragma omp parallel for                 

   for(int j = 1; j < r_b-1 ; ++j){
    //cout <<j<<endl;
    for (int i = 0; i < r_sqrt_n2; ++i){
      int index = j*r_sqrt_n2+i;
      int index2 = j*r_sqrt_n2;
      while(index2 >= M_n/2){
        index2 -= M_n/2;
      }
      batch_encoder.encode(Rot_Vec(M2[index],index2), out[index]);
      evaluator.transform_to_ntt_inplace(out[index], RLWE_ct.parms_id());
    }
  }
 // cout <<"complete 1 - "<<r_b-2<<endl;

  for (int i = 0; i < r_sqrt_n2; ++i){
      int index = (r_b-1)*r_sqrt_n2+i;
      if(index >= M_l){
        break;
      }
      int index2 = (r_b-1)*r_sqrt_n2;
      while(index2 >= M_n/2){
        index2 -= M_n/2;
      }
      batch_encoder.encode(Rot_Vec(M2[index],index2), out[index]);
      evaluator.transform_to_ntt_inplace(out[index], RLWE_ct.parms_id());
  }
 // cout <<"complete last part. "<<endl;

  
  return out;

}

extern Ciphertext LT_ecd_M_ntt(const Ciphertext & RLWE_ct, const vector<Plaintext> & M, int N,
const GaloisKeys & RotK, const SEALContext& seal_context, const SecretKey & sk){

  Evaluator evaluator(seal_context);
  

  //for test 
  Decryptor decryptor(seal_context, sk);
  BatchEncoder batch_encoder(seal_context);

  int n2 = M.size();
  size_t slot_count = batch_encoder.slot_count();
 // cout <<"size of M = "<<n2<<", slot count = "<<slot_count<<endl;

  double sqrt_n2 = sqrt((double)n2);
 // cout <<"sqrt n2 = "<<sqrt_n2<<endl;

  //rounding sqrt_n2
  int r_sqrt_n2 = (int)sqrt_n2;
  if(sqrt_n2-(double)r_sqrt_n2 > 0){
    r_sqrt_n2 ++;
  }


  //baby step

  vector<Ciphertext> c_g(r_sqrt_n2, RLWE_ct);
  evaluator.transform_to_ntt_inplace(c_g[0]);

  #pragma omp parallel for               

  for(int i = 1 ; i < r_sqrt_n2 ; ++i){
    evaluator.rotate_rows_inplace(c_g[i], i, RotK);
    evaluator.transform_to_ntt_inplace(c_g[i]);
  }

 // cout <<"Baby step end. "<<endl;

/*
  //test correctness of rotation
  cout <<"Decryption then decoding the rotated ct: "<<endl;
  for (int i = 0; i < r_sqrt_n2; ++i){
    cout <<"i = "<<i<<": ";
    Plaintext plain_result;
    decryptor.decrypt(c_g[i],plain_result);
    vector<uint64_t> pod_matrix;
    batch_encoder.decode(plain_result,pod_matrix);
    for(int j = 0 ; j < 10 ; ++j){
      cout <<pod_matrix[j]<<" ";
    }
    cout <<"...... ";
    for (int j = slot_count-10; j < slot_count; ++j){
      cout <<pod_matrix[j]<<" ";
    }
    cout <<endl;
  }
  */

  //gaint step
  double bb = (double)n2/(double)r_sqrt_n2;
  int r_b = (int)bb;
  if(bb-(double)r_b > 0){
    r_b++;
  }


  Ciphertext out(c_g[0]);
 // evaluator.multiply_plain_inplace(out,M[0]);

  //cout <<0<<endl;
  int non_zero = -1;

  //vector<Ciphertext> tempct(r_sqrt_n2);     
  //vector<bool> iszero_vec(r_sqrt_n2,true);

  for (int i = 0; i < r_sqrt_n2; ++i){
    //cout <<i<<" ";
    //check whether M[i] is zero vector
    bool iszero = true;
    for(int k = 0 ; k < slot_count ; ++k){
      if(M[i][k] != 0){
        iszero = false;
        break;
      }
    }
    if(iszero == true){
      //cout <<"M_"<<i<<" is zero vector"<<endl;
      continue;
    }
    else{
      Ciphertext tempct=c_g[i];      //tempct[i] is in ntt form
      evaluator.multiply_plain_inplace(tempct,M[i]);

      if(non_zero == -1){
           out = tempct;
           non_zero = 1;
      }
      else{
        evaluator.add_inplace(out,tempct);
      }
    }
  }

  evaluator.transform_from_ntt_inplace(out);

  vector<Ciphertext> tempout(r_b-1,c_g[0]);
  vector<int> zero_number(r_b-1,0);

  #pragma omp parallel for

  for(int j = 1; j <r_b ; ++j){
   // cout <<j<<endl;
   // int zero_number = 0;
    int first_non_zero = -1;

    //Ciphertext tempout(c_g[0]);
    //evaluator.multiply_plain_inplace(tempout,M[j*r_sqrt_n2]);

    for(int i = 0; i < r_sqrt_n2 ; ++i){
     // cout <<i<<" ";
      int index = j*r_sqrt_n2+i;
      if(index >= n2){
        break;
      }

      bool iszero = true;
      for(int k = 0 ; k < slot_count ; ++k){
        if(M[j*r_sqrt_n2+i][k] != 0){
          iszero = false;
          break;
        }
      }
      if(iszero == true){
        //cout <<"M_"<<j*r_sqrt_n2+i<<" is zero vector"<<endl;
        zero_number[j-1] ++;
        continue;
      }
      else{
        Ciphertext tempct = c_g[i];
        evaluator.multiply_plain_inplace(tempct,M[j*r_sqrt_n2+i]);
        if(first_non_zero == -1){
           tempout[j-1] = tempct;
           first_non_zero = 1;
        }
        else{
          evaluator.add_inplace(tempout[j-1],tempct);
        }
      }
      
    }

    evaluator.transform_from_ntt_inplace(tempout[j-1]);
   // cout <<endl;
    if(zero_number[j-1] == r_sqrt_n2){
      continue;
    }
    else{
      //cout <<"zero_number[j-1] != r_sqrt_n2"<<endl;
      int rot_index = j*r_sqrt_n2;
      while(rot_index >= n2/2){
        rot_index -= n2/2;
      }
      evaluator.rotate_rows_inplace(tempout[j-1],rot_index,RotK);
      
    }
    
  }

  for(int i = 1; i < r_b ; ++i){
    if(zero_number[i-1] == r_sqrt_n2){
      continue;
    }
    else{
      evaluator.add_inplace(out,tempout[i-1]);
    }
  }

  return out;

}