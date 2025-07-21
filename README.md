# HE-based-non-polynomial-eval

### Update in July
* For comparison, we include slot-encoded CKKS matrix multiplication algorithm in PCMM_CKKS folder and its test file in test folder. 

### Dependencies
* Mircosoft SEAL (attached in our thridparty folder)
* OpenMP (for multi-thread test)
* OpenSSL (for secure random number generation)
* OpenBLAS

### Installing
* Step 1: Install SEAL from thirdparty folder. Please follow the instruction of SEAL. 
* Step 2: Instal and compile
  * cd build
  * cmake ..
  * make
* Step 3: Test. There are 3 tests included. 
  * ./test
  
