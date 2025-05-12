#include "seal/seal.h"
#include "seal/util/iterator.h"
#include <numeric>
#include <stdio.h>

using namespace std;
using namespace seal;

typedef vector<int64_t> polynomial;

const int num_X = 8;
const int num_row = 128;
const int num_col = 768;
const int num_inter = 3072;
const double sqrt_d = 8.0;
int num_input = 11;
const double s0 = 16;    // s0*input
const double s1 = 1;   //s1*matrix

vector<vector<vector<double>>> input_x(num_X,vector<vector<double>>(num_row, vector<double>(num_inter,0)));
//vector<vector<double>> inter_weight(num_col, vector<double>(num_inter,0.0));
//vector<double> inter_bias(num_inter,0.0);

vector<double> gelu_plain(vector<double> &input) {
  vector<double> output;
  output.reserve(input.size());

  for (double x : input) {
    double gelu_x = 0.5 * x * (1.0 + std::tanh(std::sqrt(2.0 / M_PI) * (x + 0.044715 * x * x * x)));
    output.push_back(gelu_x);
  }

  return output;
}

void read_input(){
    ifstream fin;
    fin.open("inter_linear_output.txt");
    if(!fin.is_open()){
        cout <<"Cannot open file inter_linear_output.txt"<<endl;
    }
    char a;
    //the test file has 11 input vectors, length of each vector = 768
    
    for (int i = 0; i < num_input; ++i){
        for (int j = 0 ; j < num_inter-1 ; ++j){
            fin >>input_x[0][i][j];
            //fin >>a;
        }
        fin >>input_x[0][i][num_inter-1];
    }
    fin.close();

}


void bert_test(){
    read_input();

    omp_set_num_threads(56);
    
    cout <<"Task: gelu in bert. Input: 3072 RLWE-based ciphertexts. "<<endl;
    //Task: compute Gelu(U^T*X^T), X: 128*768, U: 768*3072
    //parameters of each group of input, X=128*768, X^T=768*128
    int input_rows = 3072;
    int input_cols = 128;
    
    cout <<"input size = "<<input_rows<<" * "<<input_cols<<endl;

    //prepare BFV
    EncryptionParameters parms(scheme_type::bfv);
    size_t poly_modulus_degree = 32768;
    parms.set_poly_modulus_degree(poly_modulus_degree);

    auto bfv_coeff_modulus = CoeffModulus::Create(poly_modulus_degree, { 30, 60, 28,55, 55, 55,
                                                          55, 55, 55, 55, 55,
                                                          55, 56 });
   // auto bfv_coeff_modulus = CoeffModulus::Create(poly_modulus_degree, { 60, 55, 28,60, 60,
   //                                                       60, 60, 60, 60, 60,
   //                                                       50, 60 });


    parms.set_coeff_modulus(bfv_coeff_modulus);
    parms.set_plain_modulus(65537);

    SEALContext context(parms);
    print_parameters(context);
    cout <<endl;
    //cout << "primitive root: " << context.first_context_data()->plain_ntt_tables()->get_root() << endl;

    auto &context_data = *context.key_context_data();
    auto coeff_modulus = context_data.parms().coeff_modulus();
    size_t coeff_modulus_size = coeff_modulus.size();
    uint64_t rlwe_ct_modulus = coeff_modulus[0].value();


    //prepare RLWE
    int RLWE_len = 1024;
    int64_t RLWE_pt_modulus = 512;
    int64_t RLWE_ct_modulus = rlwe_ct_modulus;//469762049;//576460752154525697;
    double RLWE_alpha = (double)RLWE_ct_modulus/(double)RLWE_pt_modulus;
    
    int RLWE_var = 3;

    int64_t RLWE_ct_modulus_2 = 65537;
    int64_t RLWE_ct_pt = RLWE_ct_modulus_2/RLWE_pt_modulus;

    cout <<"RLWE Parameters: Ciphertext Modulus = "<<RLWE_ct_modulus
        <<", Plaintext Modulus = "<<RLWE_pt_modulus<<", Delta in RLWE = "<<RLWE_alpha
        <<", length = "<<RLWE_len<<endl;

    int num_batch = RLWE_len/input_cols;
    cout <<"Number of input in each batch = "<<num_batch<<endl;

    auto qualifiers = context.first_context_data()->qualifiers();
    //cout << "Batching enabled: " << boolalpha << qualifiers.using_batching << endl;


    //read f from txt file
    vector<uint64_t> f(65536,1);

    
    ifstream infile;
    infile.open("lut-16.txt");
    if(!infile.is_open()){
        cout <<"cannot open file lut-16.txt. "<<endl;
    }
    int zero;
    infile >>zero;
    for (int i = 0; i < 65536; i++){
        infile >> f[i];
    }
    infile.close();
    

   // cout <<"Read f from txt file. "<<endl; 

    vector<double> div_vec(100,0.0);

    for(int iter = 0 ; iter < 1 ; ++iter){
        cout <<"-------------Round "<<iter+1<<"---------------------" <<endl;


        struct timeval tstart1, tend1;

        //prepare RLWE ct

        KeyGenerator keygen(context);
        SecretKey secret_key = keygen.secret_key();
        PublicKey public_key;
        keygen.create_public_key(public_key);
        RelinKeys relin_keys;
        keygen.create_relin_keys(relin_keys);
        Encryptor encryptor(context, public_key);
        Evaluator evaluator(context);
        Decryptor decryptor(context, secret_key);

        BatchEncoder batch_encoder(context);
        size_t slot_count = batch_encoder.slot_count();
        //cout <<"Slot count in LT = "<<slot_count<<endl;

        keygen.create_relin_keys(relin_keys);

        GaloisKeys gal_keys;
        vector<int> rot_steps = {1};
        for (int i = 0; i < RLWE_len;) {
            rot_steps.push_back(i);
            i += sqrt(RLWE_len);
        }
        keygen.create_galois_keys(rot_steps, gal_keys);

        GaloisKeys galois_keys;
        keygen.create_galois_keys(galois_keys);

        //Generate new bfv sk
        KeyGenerator new_key_keygen(context, poly_modulus_degree/RLWE_len);
        SecretKey new_key = new_key_keygen.secret_key();
        Decryptor new_decryptor(context, new_key);
        seal::util::inverse_ntt_negacyclic_harvey(new_key.data().data(), context.key_context_data()->small_ntt_tables()[0]); 

        MemoryPoolHandle my_pool = MemoryPoolHandle::New();

        // generate a key switching key based on key_before and secret_key
        KSwitchKeys ksk;
        seal::util::ConstPolyIter secret_key_before(secret_key.data().data(), poly_modulus_degree, bfv_coeff_modulus.size());

        new_key_keygen.generate_kswitch_keys(secret_key_before, 1, static_cast<KSwitchKeys &>(ksk), false);
        ksk.parms_id() = context.key_parms_id();

        // Generating RLWE secret key
        polynomial rlwe_sk(RLWE_len,0);

        int non_zero_count_after_LWE_len = 0;

        for (int i = 0; i < poly_modulus_degree; i++) {
            if(i % (poly_modulus_degree/RLWE_len) != 0){
                if(new_key.data()[i] != 0){
                    non_zero_count_after_LWE_len++;
                }
            }
            else if(new_key.data()[i] == 1){
                rlwe_sk[i/(poly_modulus_degree/RLWE_len)] = 1;
            }
            else if(new_key.data()[i] == 0){
                rlwe_sk[i/(poly_modulus_degree/RLWE_len)] = 0;
            }
            else{
                rlwe_sk[i/(poly_modulus_degree/RLWE_len)] = -1;
            }

        }

        seal::util::RNSIter new_key_rns(new_key.data().data(), poly_modulus_degree);
        ntt_negacyclic_harvey(new_key_rns, bfv_coeff_modulus.size(), context.key_context_data()->small_ntt_tables());

        vector<vector<polynomial>> output_ct(input_rows);

        for (int i = 0; i < num_inter; ++i){

            vector<double> vec(RLWE_len,0);
            for (int j = 0 ; j < 8 ; ++j){
                for (int k = 0 ; k < num_row ; ++k){
                    vec[j*num_row+k] = s0*input_x[0][k][i];
                }
            }

            output_ct[i] = RLWE64_Enc(RLWE_len, RLWE_ct_modulus, RLWE_alpha, RLWE_var, vec, rlwe_sk);
        }

        //decryption test
        cout <<"FOR TEST: decryption result of first 10 RLWE cts: "<<endl;
        for (int i = 0; i < 10; ++i){
            vector<double> dec_m = RLWE64_Dec(RLWE_len, RLWE_ct_modulus, RLWE_alpha, rlwe_sk, output_ct[i]);
            for (int j = 0 ; j < 11 ; ++j){
            cout <<dec_m[j]/s0<<" ";
            }
            cout<<endl;
        }
        cout <<endl;

        cout <<"======= TEST START ========"<<endl;

        //Generating extracted ct
        int output_num = output_ct.size();
        cout <<"number of input ct = "<<output_num<<endl;
        int extracted_ct_num = output_num*input_cols*num_batch;
        cout <<"total number of message = "<<extracted_ct_num<<endl;

        gettimeofday(&tstart1,NULL);

        for (int i = 0; i < output_num; ++i){
            output_ct[i] = RLWE64_Rounding(RLWE_ct_modulus,RLWE_ct_modulus_2,RLWE_len,output_ct[i]);
        }

        vector<vector<uint64_t>> extracted_ct = Extract_rlwe(output_ct,output_num,RLWE_len,RLWE_len,RLWE_ct_modulus_2);

        vector<vector<uint64_t>> RLWE_ct2 = reorder(extracted_ct,output_num,RLWE_len,RLWE_len);

        double RLWE_alpha_2 = RLWE_alpha/(double)(RLWE_ct_modulus/RLWE_ct_modulus_2);

        gettimeofday(&tend1,NULL);
        double ext_time = tend1.tv_sec-tstart1.tv_sec+(tend1.tv_usec-tstart1.tv_usec)/1000000.0;

        //LT to BFV
        //prepare RLWE(ecd(LWE_sk))
        Ciphertext rlwe_sk_encrypted = encryptLWEskUnderBFV(context, poly_modulus_degree,
            public_key, secret_key, rlwe_sk, RLWE_len, RLWE_ct_modulus_2);

        //construct matrix M
        gettimeofday(&tstart1,NULL);

        int sq_sk = sqrt(RLWE_len);
        vector<Ciphertext> rlwe_sk_sqrt_list(sq_sk);

        Ciphertext rlwe_sk_column;

        auto context_data_ptr = context.get_context_data(rlwe_sk_encrypted.parms_id());

        evaluator.rotate_columns(rlwe_sk_encrypted, gal_keys, rlwe_sk_column);

        #pragma omp parallel for           
        for (int i = 0; i < sq_sk; i++) {
            evaluator.rotate_rows(rlwe_sk_encrypted, sq_sk * i, gal_keys, rlwe_sk_sqrt_list[i]);
            evaluator.transform_to_ntt_inplace(rlwe_sk_sqrt_list[i]);
        }

        vector<Ciphertext> encrypted_result = evaluatePackedLWECiphertext_batch(context, RLWE_ct2, rlwe_sk_sqrt_list, 
            gal_keys, RLWE_len, poly_modulus_degree, RLWE_ct_modulus_2);

        gettimeofday(&tend1,NULL);
        double LT_time = tend1.tv_sec-tstart1.tv_sec+(tend1.tv_usec-tstart1.tv_usec)/1000000.0;

        int num_bfv_ct = encrypted_result.size();

        cout <<"Time for RLWE -> BFV = "<<LT_time+ext_time<<"s. "<<endl;
        cout << "Noise budget in LT result: " << decryptor.invariant_noise_budget(encrypted_result[0]) << " bits"<< endl;
        cout <<"number of primes in coeff_modulus of LT result: "<<encrypted_result[0].coeff_modulus_size()<<endl;
        cout <<"Number of BFV ciphertexts = "<<num_bfv_ct;

        
        gettimeofday(&tstart1,NULL);

        //RLWE(ecd(f(Ask+b)))
        vector<Ciphertext> encrypted_eval_result(num_bfv_ct);

        #pragma omp parallel for 

        for (int i = 0; i < num_bfv_ct; ++i){
            Bootstrap_RangeCheck_PatersonStockmeyer(encrypted_eval_result[i], encrypted_result[i], f, poly_modulus_degree, relin_keys, context, secret_key);
            evaluator.mod_switch_to_next_inplace(encrypted_eval_result[i]);
        }
        
        gettimeofday(&tend1,NULL);
        double poly_eval_time = tend1.tv_sec-tstart1.tv_sec+(tend1.tv_usec-tstart1.tv_usec)/1000000.0;
        
        cout <<"Time for polynomial evaluation = "<<poly_eval_time<<"s. "<<endl;
        cout <<"number of primes in coeff_modulus of polynomial evaluation: "<<encrypted_eval_result[0].coeff_modulus_size()<<endl;

        //ntt
        evaluator.transform_to_ntt_inplace(encrypted_eval_result[0]);

        gettimeofday(&tstart1,NULL);

        //construct dcd matrix
        vector<vector<uint64_t>> dcdmatrix = cal_dcd_matrix(slot_count,context);

        vector<Plaintext> ecd_dcdmatrix = ecd_dcd_matrix(poly_modulus_degree,dcdmatrix,context,encrypted_eval_result[0]);

        gettimeofday(&tend1,NULL);
        double pre_pro_dcd_matrix = tend1.tv_sec-tstart1.tv_sec+(tend1.tv_usec-tstart1.tv_usec)/1000000.0;
        // cout <<"Identity function evaluation end. ";
        cout <<"Time for pre-processing decoding matrix = "<<pre_pro_dcd_matrix<<"s. "<<endl;

        //intt
        evaluator.transform_from_ntt_inplace(encrypted_eval_result[0]);

        //cout <<"S2C start"<<endl;
        vector<Ciphertext> LT_dcd_out_1(num_bfv_ct);
  
        gettimeofday(&tstart1,NULL);

        #pragma omp parallel for

        for(int i = 0 ; i < num_bfv_ct ; ++i){
            LT_dcd_out_1[i] = Slot_to_coeff(encrypted_eval_result[i], ecd_dcdmatrix, galois_keys, context, secret_key);
        }
        
        gettimeofday(&tend1,NULL);
        double dcdLT_time = tend1.tv_sec-tstart1.tv_sec+(tend1.tv_usec-tstart1.tv_usec)/1000000.0;
        cout <<"Time for decoding LT = "<<dcdLT_time<<"s. "<<endl;
        cout <<"number of primes in coeff_modulus of decoding LT: "<<LT_dcd_out_1[0].coeff_modulus_size()<<endl;


        gettimeofday(&tstart1,NULL);

        //mod switch to smaller modulus
        #pragma omp parallel for

        for (int i = 0; i < num_bfv_ct; ++i){
            while(context.last_parms_id() != LT_dcd_out_1[i].parms_id()){
                evaluator.mod_switch_to_next_inplace(LT_dcd_out_1[i]);
            }
        }
        
        #pragma omp parallel for

        for (int i = 0; i < num_bfv_ct; ++i){
            Ciphertext copy_coeff = LT_dcd_out_1[i];
            auto ct_in_iter = util::iter(copy_coeff);
            ct_in_iter += LT_dcd_out_1[i].size() - 1;
            seal::util::set_zero_poly(poly_modulus_degree, 1, LT_dcd_out_1[i].data(1)); // notice that the coeff_mod.size() is hardcoded to 1, thus this needs to be performed on the last level
            evaluator.switch_key_inplace(LT_dcd_out_1[i], *ct_in_iter, static_cast<const KSwitchKeys &>(ksk), 0, my_pool);
        }
        

        gettimeofday(&tend1,NULL);
        double ksk_time = tend1.tv_sec-tstart1.tv_sec+(tend1.tv_usec-tstart1.tv_usec)/1000000.0;
        cout <<"Time for key switch + mod switch = "<<ksk_time<<"s. "<<endl;


        // TO DO: Ring switch to 1024-degree polynomial ring
        gettimeofday(&tstart1,NULL);

        vector<vector<polynomial>> final_rlwe_ct = Ring_switch(LT_dcd_out_1,poly_modulus_degree,RLWE_len);

        gettimeofday(&tend1,NULL);
        double rsk_time = tend1.tv_sec-tstart1.tv_sec+(tend1.tv_usec-tstart1.tv_usec)/1000000.0;
        cout <<"Time for ring switch = "<<rsk_time<<"s. "<<endl;
        cout <<"number of ciphertexts = "<<final_rlwe_ct.size();

        vector<vector<double>> result(final_rlwe_ct.size(),vector<double>(RLWE_len,0));

        //decryption test
        cout <<"FOR TEST: decryption result of first 10 RLWE cts: "<<endl;
        double error = 0;
        for (int i = 0; i < 10; ++i){
            result[i] = RLWE64_Dec(RLWE_len, rlwe_ct_modulus, (double)rlwe_ct_modulus/(double)RLWE_pt_modulus, 
                rlwe_sk, final_rlwe_ct[i]);
            //if(i <= 0){
                for (int j = 0 ; j < 11 ; ++j){
                    cout <<result[i][j]/s0<<" ";
                }
                cout <<endl;
            //}

        }


 
    }
    

}