using namespace std;
using namespace seal;

vector<Ciphertext> ct_pt_matrix_mul_wo_pre(const vector<Ciphertext> & enc_X, 
  const vector<vector<double>> & W, int col_X, int col_W, int row_W, 
  const SEALContext& seal_context){

  vector<Ciphertext> output(col_W);
  double scale = enc_X[0].scale();

  if(col_X != row_W){
    cout <<"ERROR: bad dimensions of X or W. "<<endl;
    return output;
  }

  CKKSEncoder encoder(seal_context);
  // Evaluator evaluator(seal_context, encoder);
  Evaluator evaluator(seal_context);

    #pragma omp parallel for 
    for (int i = 0; i < col_W; ++i){
        //encode w[0][i]
        Plaintext ecd_w_0_i;
        encoder.encode(W[0][i], enc_X[0].parms_id(), enc_X[0].scale(), ecd_w_0_i);
        //enc_X[0]*ecd_w[0][i]
        evaluator.multiply_plain(enc_X[0], ecd_w_0_i, output[i]);
        //evaluator.rescale_to_next_inplace(output[i]);

        for (int j = 1 ; j < row_W ; ++j){
          //encode w[j][i]
          Plaintext ecd_w_j_i;
          encoder.encode(W[j][i], enc_X[j].parms_id(), enc_X[j].scale(), ecd_w_j_i);

          //enc_X[j]*ecd_w[j][i]
          Ciphertext temp;
          evaluator.multiply_plain(enc_X[j], ecd_w_j_i, temp);
          //evaluator.rescale_to_next_inplace(temp);
          evaluator.add_inplace(output[i],temp);
        }

        evaluator.rescale_to_next_inplace(output[i]);
        output[i].scale()=scale;
    }


  

  return output;

}

