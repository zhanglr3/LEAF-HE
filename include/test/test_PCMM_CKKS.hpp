using namespace std;
using namespace seal;

void ct_pt_matrix_mul_test(){
    cout <<"Task: test slot-encoded Ct-Pt matrix multiplication in CKKS scheme: "<<endl;

    EncryptionParameters parms(scheme_type::ckks);

    size_t poly_modulus_degree = 1 << 16;
    parms.set_poly_modulus_degree(poly_modulus_degree);
   // parms.set_coeff_modulus(CoeffModulus::BFVDefault(poly_modulus_degree));
    parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, 
        {60, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 60}));
   //     {60, 40, 40, 60}));
    double scale = pow(2.0,40);

    //SEALContext context(parms);
    SEALContext context(parms, true, sec_level_type::none);

    cout <<"Set encryption parameters and print"<<endl;
    print_parameters(context);

    KeyGenerator keygen(context);
    SecretKey secret_key = keygen.secret_key();
    PublicKey public_key;
    keygen.create_public_key(public_key);

    Decryptor decryptor(context, secret_key);
    CKKSEncoder encoder(context);
    Evaluator evaluator(context);
    size_t slot_count = encoder.slot_count();

    struct timeval tstart1, tend1;

    //construct input
    int num_X = 256;
    int num_row = 128;
    int num_col = 768;
    vector<vector<vector<double>>> input_x(num_X,vector<vector<double>>(num_row, vector<double>(num_col,0)));
    for (int i = 0; i < num_X; ++i){
        for (int j = 0 ; j < num_row ; ++j){
            for (int k = 0 ; k < num_col ; ++k){
                input_x[i][j][k] = (double)j+1.0;
            }
        }

    }
    cout <<"Matrix X size = "<<num_row <<" * "<<num_col<<endl;
    

    //encode + encrypt
    vector<Ciphertext> enc_ecd_x = batch_input(input_x, num_X, num_row, num_col, scale, context,public_key);

    cout <<"encode and encrypt X. "<<endl;
    cout <<"Modulus chain index for enc x: "<< context.get_context_data(enc_ecd_x[0].parms_id())->chain_index()<<endl;

    //construct W
    int col_W = 3072;
    vector<vector<double>> W(num_col, vector<double>(col_W, 1.0/128.0));
    cout <<"Matrix W size = "<<num_col <<" * "<<col_W<<endl;

    cout <<"Encrypted col-packing X * W = Encrypted col-packing XW. "<<endl;

    //matrix multiplication
    gettimeofday(&tstart1,NULL);

    vector<Ciphertext> ct_pt_mul = ct_pt_matrix_mul_wo_pre(enc_ecd_x, W, num_col, col_W, num_col, context);

    gettimeofday(&tend1,NULL);
    double ct_pt_matrix_mul_time = tend1.tv_sec-tstart1.tv_sec+(tend1.tv_usec-tstart1.tv_usec)/1000000.0;
    cout <<"Ct-Pt matrix multiplication time = "<<ct_pt_matrix_mul_time<<endl;
    cout <<"Modulus chain index for the result: "<< context.get_context_data(ct_pt_mul[0].parms_id())->chain_index()<<endl;

    cout <<"Decrypt + decode result: "<<endl;
    //decrypt and decode
    for (int i = 0; i < 5; ++i){
        Plaintext plain_result;
        decryptor.decrypt(ct_pt_mul[i], plain_result);
        vector<double> result;
        encoder.decode(plain_result, result);
        cout <<i+1<<"-th ciphertext: ";
        for (int ind = 0 ; ind < 5 ; ++ind){
            cout <<result[ind]<<" ";
        }
        cout <<"... ";
        for (int ind = slot_count-5 ; ind < slot_count ; ++ind){
            cout <<result[ind]<<" ";
        }
        cout <<endl;
    }
    cout <<"......"<<endl;
    for (int i = col_W-5; i < col_W; ++i){
        Plaintext plain_result;
        decryptor.decrypt(ct_pt_mul[i], plain_result);
        vector<double> result;
        encoder.decode(plain_result, result);
        cout <<i+1<<"-th ciphertext: ";
        for (int ind = 0 ; ind < 5 ; ++ind){
            cout <<result[ind]<<" ";
        }
        cout <<"... ";
        for (int ind = slot_count-5 ; ind < slot_count ; ++ind){
            cout <<result[ind]<<" ";
        }
        cout <<endl;
    }
}
