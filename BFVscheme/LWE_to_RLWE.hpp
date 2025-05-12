Ciphertext encryptLWEskUnderBFV(const SEALContext& context, const size_t& degree,
                                const PublicKey& BFVpk, const SecretKey& BFVsk,
                                const vector<int64_t>& LWEsk, const int LWE_len, const int64_t LWE_ct_mod) { 
    Ciphertext switchingKey(context);

    BatchEncoder batch_encoder(context);
    Encryptor encryptor(context, BFVpk);
    encryptor.set_secret_key(BFVsk);

    int tempn = 1;
    for(tempn = 1; tempn < LWE_len; tempn *= 2){}

    vector<uint64_t> skInt(degree);
    for(size_t i = 0; i < degree; i++){
        auto tempindex = i%uint64_t(tempn);
        if(int(tempindex) >= LWE_len) {
            skInt[i] = 0;
        } else {
        	if (LWEsk[tempindex] >= 0){
        		skInt[i] = (uint64_t)LWEsk[tempindex];
        	}
        	else{
        		skInt[i] = uint64_t( LWE_ct_mod - 1);
        	}
            
        }
    }
    Plaintext plaintext;
    batch_encoder.encode(skInt, plaintext);
    encryptor.encrypt_symmetric(plaintext, switchingKey);

    return switchingKey;
}

// assume lwe_sk_len is a power of 2, and has a square root
Ciphertext evaluatePackedLWECiphertext(const SEALContext& seal_context, vector<vector<uint64_t> >& lwe_ct_list,
                                       const vector<Ciphertext>& lwe_sk_sqrt_list, const GaloisKeys& gal_keys, const int lwe_sk_len,
                                       const int degree, const int q = 65537) {
    Evaluator evaluator(seal_context);
    BatchEncoder batch_encoder(seal_context);

    // rotate sqrt(degree), and get sqrt(degree)'s lwe_sk_encrypted
    int sq_rt = sqrt(lwe_sk_len);
    vector<Ciphertext> result(sq_rt);

    #pragma omp parallel for               
        
    for (int iter = 0; iter < sq_rt; iter++) {
        for (int j = 0; j < (int) lwe_sk_sqrt_list.size(); j++) {
            vector<uint64_t> lwe_ct_tmp(degree);
            for (int i = 0; i < degree; i++) {
                int ct_index = (i-iter) % (degree/2) < 0 ? (i-iter) % (degree/2) + degree/2 : (i-iter) % (degree/2);
                ct_index = i < degree/2 ? ct_index : ct_index + degree/2;
                int col_index = (i + j*sq_rt) % lwe_sk_len;
                lwe_ct_tmp[i] = lwe_ct_list[ct_index][col_index];
            }

            Plaintext lwe_ct_pl;
            batch_encoder.encode(lwe_ct_tmp, lwe_ct_pl);
            evaluator.transform_to_ntt_inplace(lwe_ct_pl, lwe_sk_sqrt_list[j].parms_id());

            if (j == 0) {
                evaluator.multiply_plain(lwe_sk_sqrt_list[j], lwe_ct_pl, result[iter]);
            } else {
                Ciphertext temp;
                evaluator.multiply_plain(lwe_sk_sqrt_list[j], lwe_ct_pl, temp);
                evaluator.add_inplace(result[iter], temp);
            }

        }
    }

    for (int i = 0; i < (int) result.size(); i++) {
        evaluator.transform_from_ntt_inplace(result[i]);
    }

    // sum up all sq_rt tmp results to the first one, each first rotate left one and add to the previous tmp result

    for (int iter = sq_rt-1; iter > 0; iter--) {
        evaluator.rotate_rows_inplace(result[iter], 1, gal_keys);
        evaluator.add_inplace(result[iter-1], result[iter]);
    }

    vector<uint64_t> b_parts(degree);
    for(int i = 0; i < degree; i++){
        b_parts[i] = lwe_ct_list[i][lwe_sk_len];
    }

    Plaintext lwe_b_pl;
    batch_encoder.encode(b_parts, lwe_b_pl);
    //evaluator.negate_inplace(result[0]);
    evaluator.add_plain_inplace(result[0], lwe_b_pl);

    return result[0];
}

// assume lwe_sk_len is a power of 2, and has a square root
vector<Ciphertext> evaluatePackedLWECiphertext_batch(const SEALContext& seal_context, vector<vector<uint64_t> >& lwe_ct_list,
                                       const vector<Ciphertext>& lwe_sk_sqrt_list, const GaloisKeys& gal_keys, const int lwe_sk_len,
                                       const int degree, const int q = 65537) {
    Evaluator evaluator(seal_context);
    BatchEncoder batch_encoder(seal_context);

    // rotate sqrt(degree), and get sqrt(degree)'s lwe_sk_encrypted
    int sq_rt = sqrt(lwe_sk_len);

    int num_bfv_ct = lwe_ct_list.size()/degree;
    if(num_bfv_ct * degree < lwe_ct_list.size()){
        num_bfv_ct ++;
    }
    //cout <<"number of BFV ct: "<<num_bfv_ct<<endl;
    vector<Ciphertext> output(num_bfv_ct);

    #pragma omp parallel for 

    for(int ind = 0 ; ind < num_bfv_ct ; ++ind){

        vector<Ciphertext> result(sq_rt);

        for (int iter = 0; iter < sq_rt; iter++) {
            for (int j = 0; j < (int) lwe_sk_sqrt_list.size(); j++) {
                vector<uint64_t> lwe_ct_tmp(degree);
                for (int i = 0; i < degree; i++) {
                    int ct_index = (i-iter) % (degree/2) < 0 ? (i-iter) % (degree/2) + degree/2 : (i-iter) % (degree/2);
                    ct_index = i < degree/2 ? ct_index : ct_index + degree/2;
                    int col_index = (i + j*sq_rt) % lwe_sk_len;
                    lwe_ct_tmp[i] = lwe_ct_list[ct_index+ind*degree][col_index];
                }

                Plaintext lwe_ct_pl;
                batch_encoder.encode(lwe_ct_tmp, lwe_ct_pl);
                evaluator.transform_to_ntt_inplace(lwe_ct_pl, lwe_sk_sqrt_list[j].parms_id());

                if (j == 0) {
                    evaluator.multiply_plain(lwe_sk_sqrt_list[j], lwe_ct_pl, result[iter]);
                } 
                else {
                    Ciphertext temp;
                    evaluator.multiply_plain(lwe_sk_sqrt_list[j], lwe_ct_pl, temp);
                    evaluator.add_inplace(result[iter], temp);
                }

            }
        }

        for (int i = 0; i < (int) result.size(); i++) {
            evaluator.transform_from_ntt_inplace(result[i]);
        }

        // sum up all sq_rt tmp results to the first one, each first rotate left one and add to the previous tmp result

        for (int iter = sq_rt-1; iter > 0; iter--) {
            evaluator.rotate_rows_inplace(result[iter], 1, gal_keys);
            evaluator.add_inplace(result[iter-1], result[iter]);
        }

        vector<uint64_t> b_parts(degree);
        for(int i = 0; i < degree; i++){
            b_parts[i] = lwe_ct_list[i+ind*degree][lwe_sk_len];
        }

        Plaintext lwe_b_pl;
        batch_encoder.encode(b_parts, lwe_b_pl);
        //evaluator.negate_inplace(result[0]);
        evaluator.add_plain_inplace(result[0], lwe_b_pl);
        output[ind] = result[0];
    }

    return output;
}

