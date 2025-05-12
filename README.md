# HE-based-non-polynomial-eval

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
  
