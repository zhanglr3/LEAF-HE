#include "seal/seal.h"
#include "seal/util/iterator.h"
#include <numeric>
#include <stdio.h>

using namespace std;
using namespace seal;

typedef vector<int64_t> polynomial;

inline void print_parameters(const seal::SEALContext &context)
{
    auto &context_data = *context.key_context_data();

    /*
    Which scheme are we using?
    */
    std::string scheme_name;
    switch (context_data.parms().scheme())
    {
    case seal::scheme_type::bfv:
        scheme_name = "BFV";
        break;
    case seal::scheme_type::ckks:
        scheme_name = "CKKS";
        break;
    case seal::scheme_type::bgv:
        scheme_name = "BGV";
        break;
    default:
        throw std::invalid_argument("unsupported scheme");
    }
    std::cout << "/" << std::endl;
    std::cout << "| Encryption parameters :" << std::endl;
    std::cout << "|   scheme: " << scheme_name << std::endl;
    std::cout << "|   poly_modulus_degree: " << context_data.parms().poly_modulus_degree() << std::endl;

    /*
    Print the size of the true (product) coefficient modulus.
    */
    std::cout << "|   coeff_modulus size: ";
    std::cout << context_data.total_coeff_modulus_bit_count() << " (";
    auto coeff_modulus = context_data.parms().coeff_modulus();
    std::size_t coeff_modulus_size = coeff_modulus.size();
    for (std::size_t i = 0; i < coeff_modulus_size - 1; i++)
    {
        std::cout << coeff_modulus[i].bit_count() << " + ";
    }
    std::cout << coeff_modulus.back().bit_count();
    std::cout << ") bits" << std::endl;

    /*
    For the BFV scheme print the plain_modulus parameter.
    */
    if (context_data.parms().scheme() == seal::scheme_type::bfv)
    {
        std::cout << "|   plain_modulus: " << context_data.parms().plain_modulus().value() << std::endl;
    }

    std::cout << "\\" << std::endl;
}


void full_protocol_test(){
    cout <<"Task: Full protocol. "<<endl;
    //parameters of each group of input
    int input_rows = 128;
    int input_cols = 1024;
    int input_batch = 1;
    cout <<"input size = "<<input_rows<<" * "<<input_cols<<", number of batch = "<<input_batch<<endl;

    //parameters of matrix U
    int matrix_rows = 32;
    int matrix_cols = input_rows;
    cout <<"matrix size = "<<matrix_rows <<" * "<<matrix_cols<<endl;


    

    //prepare BFV
    EncryptionParameters parms(scheme_type::bfv);
    size_t poly_modulus_degree = 32768;
    parms.set_poly_modulus_degree(poly_modulus_degree);

    auto bfv_coeff_modulus = CoeffModulus::Create(poly_modulus_degree, { 30, 60, 28,55, 55, 55,
                                                          55, 55, 55, 55, 55,
                                                          55, 56 });

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
    //cout <<rlwe_ct_modulus<<endl;

    //prepare RLWE
    int RLWE_len = 1024;
    int64_t RLWE_pt_modulus = 512;
    int64_t RLWE_ct_modulus = rlwe_ct_modulus;
    double RLWE_alpha = (double)RLWE_ct_modulus/(double)RLWE_pt_modulus;

    int64_t RLWE_ct_modulus_2 = 65537;
    int64_t RLWE_ct_pt = RLWE_ct_modulus_2/RLWE_pt_modulus;

    
    int RLWE_var = 3;

    cout <<"RLWE Parameters: Ciphertext Modulus = "<<RLWE_ct_modulus
        <<", Plaintext Modulus = "<<RLWE_pt_modulus<<", Delta in RLWE = "<<RLWE_alpha
        <<", length = "<<RLWE_len<<endl;

/*
    //prepare LWE
    int64_t LWE_pt_modulus = 512;
    int64_t LWE_ct_modulus = 65537;
    int64_t LWE_ct_pt = LWE_ct_modulus/LWE_pt_modulus;
    int64_t ecd_modulus = LWE_ct_modulus;
    int LWE_len = 1024;
    int LWE_var = 3;
*/
    auto qualifiers = context.first_context_data()->qualifiers();
    //cout << "Batching enabled: " << boolalpha << qualifiers.using_batching << endl;


    //read f from txt file
    vector<uint64_t> f(65536,1);

    
    ifstream infile;
    infile.open("fpoly-65537-512.txt");
    if(!infile.is_open()){
        cout <<"cannot open file fpoly-65537-512.txt. "<<endl;
    }
    int zero;
    infile >>zero;
    for (int i = 0; i < 65536; i++){
        infile >> f[i];
    }
    infile.close();
    

   // cout <<"Read f from txt file. "<<endl; 

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
                //cout << i/(poly_modulus_degree/RLWE_len)<<" ";
                rlwe_sk[i/(poly_modulus_degree/RLWE_len)] = 1;
            }
            else if(new_key.data()[i] == 0){
                //cout << i/(poly_modulus_degree/RLWE_len)<<" ";
                rlwe_sk[i/(poly_modulus_degree/RLWE_len)] = 0;
            }
            else{
                //cout << i/(poly_modulus_degree/RLWE_len)<<" ";
                rlwe_sk[i/(poly_modulus_degree/RLWE_len)] = -1;
            }

        }

        seal::util::RNSIter new_key_rns(new_key.data().data(), poly_modulus_degree);
        ntt_negacyclic_harvey(new_key_rns, bfv_coeff_modulus.size(), context.key_context_data()->small_ntt_tables());

        //Generate RLWE ct
        vector<double> m(RLWE_len,0);
        for (int i = 0; i < input_cols; ++i){
            m[i] = (double)(i+1)*0.01;
        }

        vector<vector<polynomial>> input_ct(input_rows);
        for (int i = 0; i < input_rows; ++i){
            input_ct[i] = RLWE64_Enc(RLWE_len, RLWE_ct_modulus, RLWE_alpha, RLWE_var, m, rlwe_sk);
        }

        //Generate matrix U
        vector<vector<float>> matrixU(matrix_rows,vector<float>(matrix_cols,1));
        for (int i = 0; i < matrix_rows; ++i){
            for (int j = 0 ; j < matrix_cols  ; ++j){
                matrixU[i][j] = j % 10;
            }
        }

        //prepare RLWE(ecd(LWE_sk))
        Ciphertext rlwe_sk_encrypted = encryptLWEskUnderBFV(context, poly_modulus_degree,
            public_key, secret_key, rlwe_sk, RLWE_len, RLWE_ct_modulus_2);

        cout <<"======= TEST START ========"<<endl;
        //Compute U*RLWE ct
        gettimeofday(&tstart1,NULL);
        vector<vector<polynomial>> output_ct = RLWE64_matrix_mul(matrix_rows,matrix_cols,RLWE_len,RLWE_ct_modulus,input_ct,matrixU);
        int output_num = output_ct.size();
        gettimeofday(&tend1,NULL);
        double mm_time = tend1.tv_sec-tstart1.tv_sec+(tend1.tv_usec-tstart1.tv_usec)/1000000.0;
        cout <<"Matrix multiplication. number of output ct = "<<output_ct.size()<<endl;
        cout <<"Time for matrix multiplication = "<<mm_time<<"s. "<<endl;

        //decryption test
        cout <<"FOR TEST: decryption result of first 10 RLWE cts: "<<endl;
        for (int i = 0; i < 10; ++i){
            vector<double> dec_m = RLWE64_Dec(RLWE_len, RLWE_ct_modulus, RLWE_alpha, rlwe_sk, output_ct[i]);
            for (int j = 0 ; j < 10 ; ++j){
            cout <<dec_m[j]<<" ";
            }
            cout <<endl;
        }
        cout <<endl;


        //Generating LWE ct
        int LWE_ct_num = output_num*input_cols;
        //cout <<"total number of message = "<<LWE_ct_num<<endl;

        gettimeofday(&tstart1,NULL);
        for (int i = 0; i < output_num; ++i){
            output_ct[i] = RLWE64_Rounding(RLWE_ct_modulus,RLWE_ct_modulus_2,RLWE_len,output_ct[i]);
        }
        vector<vector<uint64_t>> extracted_ct = Extract_rlwe(output_ct,output_num,RLWE_len,RLWE_len,RLWE_ct_modulus_2);

        vector<vector<uint64_t>> RLWE_ct2 = reorder_single(extracted_ct,output_num,RLWE_len,RLWE_len);

        double RLWE_alpha_2 = RLWE_alpha/(double)(RLWE_ct_modulus/RLWE_ct_modulus_2);

        //LT to BFV

        int sq_sk = sqrt(RLWE_len);
        vector<Ciphertext> rlwe_sk_sqrt_list(sq_sk);
        //cout <<"sq_sk = "<<sq_sk<<endl;

        Ciphertext rlwe_sk_column;

        auto context_data_ptr = context.get_context_data(rlwe_sk_encrypted.parms_id());

        evaluator.rotate_columns(rlwe_sk_encrypted, gal_keys, rlwe_sk_column);

        #pragma omp parallel for           
        for (int i = 0; i < sq_sk; i++) {
            evaluator.rotate_rows(rlwe_sk_encrypted, sq_sk * i, gal_keys, rlwe_sk_sqrt_list[i]);
            evaluator.transform_to_ntt_inplace(rlwe_sk_sqrt_list[i]);
        }

        Ciphertext encrypted_result = evaluatePackedLWECiphertext(context, RLWE_ct2, rlwe_sk_sqrt_list, 
            gal_keys, RLWE_len, poly_modulus_degree, RLWE_ct_modulus_2);

        gettimeofday(&tend1,NULL);
        double LT_time = tend1.tv_sec-tstart1.tv_sec+(tend1.tv_usec-tstart1.tv_usec)/1000000.0;

        cout <<"Time for RLWE -> BFV = "<<LT_time<<"s. "<<endl;
        cout << "Noise budget in LT result: " << decryptor.invariant_noise_budget(encrypted_result) << " bits"<< endl;
        cout <<"number of primes in coeff_modulus of LT result: "<<encrypted_result.coeff_modulus_size()<<endl;
        cout <<endl;

        //test correctness of LT
        Plaintext plain_result;
        decryptor.decrypt(encrypted_result, plain_result);
        vector<uint64_t> dcd_matrix(slot_count,0);
        batch_encoder.decode(plain_result, dcd_matrix);


        gettimeofday(&tstart1,NULL);

        //RLWE(ecd(f(Ask+b)))
        Ciphertext encrypted_eval_result;
        Bootstrap_RangeCheck_PatersonStockmeyer(encrypted_eval_result, encrypted_result, f, poly_modulus_degree, relin_keys, context, secret_key);
        evaluator.mod_switch_to_next_inplace(encrypted_eval_result);
        gettimeofday(&tend1,NULL);
        double poly_eval_time = tend1.tv_sec-tstart1.tv_sec+(tend1.tv_usec-tstart1.tv_usec)/1000000.0;
        
        cout <<"Time for polynomial evaluation = "<<poly_eval_time<<"s. "<<endl;
        cout <<"number of primes in coeff_modulus of polynomial evaluation: "<<encrypted_eval_result.coeff_modulus_size()<<endl;

        Plaintext plain_result2;
        decryptor.decrypt(encrypted_eval_result, plain_result2);
        vector<uint64_t> dcd_matrix2(slot_count,0);
        batch_encoder.decode(plain_result2, dcd_matrix2);


        //ntt
        evaluator.transform_to_ntt_inplace(encrypted_eval_result);

        gettimeofday(&tstart1,NULL);

        //construct dcd matrix
        vector<vector<uint64_t>> dcdmatrix = cal_dcd_matrix(slot_count,context);

        vector<Plaintext> ecd_dcdmatrix = ecd_dcd_matrix(poly_modulus_degree,dcdmatrix,context,encrypted_eval_result);


        //cout <<"decoding matrix encoded. "<<endl;

        gettimeofday(&tend1,NULL);
        double pre_pro_dcd_matrix = tend1.tv_sec-tstart1.tv_sec+(tend1.tv_usec-tstart1.tv_usec)/1000000.0;
        // cout <<"Identity function evaluation end. ";
        cout <<"Time for pre-processing decoding matrix = "<<pre_pro_dcd_matrix<<"s. "<<endl;

        //intt
        evaluator.transform_from_ntt_inplace(encrypted_eval_result);

        //cout <<"S2C start"<<endl;
  
        gettimeofday(&tstart1,NULL);

        Ciphertext LT_dcd_out_1 = Slot_to_coeff(encrypted_eval_result, ecd_dcdmatrix, galois_keys, context, secret_key);

        gettimeofday(&tend1,NULL);
        double dcdLT_time = tend1.tv_sec-tstart1.tv_sec+(tend1.tv_usec-tstart1.tv_usec)/1000000.0;
        cout <<"Time for decoding LT = "<<dcdLT_time<<"s. "<<endl;
        cout <<"number of primes in coeff_modulus of decoding LT: "<<LT_dcd_out_1.coeff_modulus_size()<<endl;

        Plaintext plain_result3;
        decryptor.decrypt(LT_dcd_out_1, plain_result3);


        gettimeofday(&tstart1,NULL);

        //mod switch to smaller modulus
        while(context.last_parms_id() != LT_dcd_out_1.parms_id()){
            evaluator.mod_switch_to_next_inplace(LT_dcd_out_1);
        }

        
        Ciphertext copy_coeff = LT_dcd_out_1;
        auto ct_in_iter = util::iter(copy_coeff);
        ct_in_iter += LT_dcd_out_1.size() - 1;
        seal::util::set_zero_poly(poly_modulus_degree, 1, LT_dcd_out_1.data(1)); // notice that the coeff_mod.size() is hardcoded to 1, thus this needs to be performed on the last level

        evaluator.switch_key_inplace(LT_dcd_out_1, *ct_in_iter, static_cast<const KSwitchKeys &>(ksk), 0, my_pool);

        gettimeofday(&tend1,NULL);
        double ksk_time = tend1.tv_sec-tstart1.tv_sec+(tend1.tv_usec-tstart1.tv_usec)/1000000.0;
        cout <<"Time for key switch + mod switch = "<<ksk_time<<"s. "<<endl;


        // TO DO: Ring switch to 1024-degree polynomial ring
        gettimeofday(&tstart1,NULL);

        vector<vector<polynomial>> final_rlwe_ct = Ring_switch_single(LT_dcd_out_1,poly_modulus_degree,RLWE_len);

        gettimeofday(&tend1,NULL);
        double ring_switch_time = tend1.tv_sec-tstart1.tv_sec+(tend1.tv_usec-tstart1.tv_usec)/1000000.0;
        cout <<"Time for ring switch = "<<ring_switch_time<<"s. "<<endl;

        //decryption test
        cout <<"FOR TEST: decryption result of first 10 RLWE cts: "<<endl;
        for (int i = 0; i < 10; ++i){
            vector<double> dec_m = RLWE64_Dec(RLWE_len, rlwe_ct_modulus, (double)rlwe_ct_modulus/(double)RLWE_pt_modulus, 
                rlwe_sk, final_rlwe_ct[i]);
            for (int j = 0 ; j < 10 ; ++j){
            cout <<dec_m[j]<<" ";
            }
            cout <<endl;
        }
        cout <<endl;

 
    }
    

}