#ifndef INCLUDE_H
#define INCLUDE_H

#pragma once

//SEAL library
#include "seal/seal.h"

//C++
#include <iostream>
#include <vector>
#include <ctime>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>
#include <fstream>
#include <cstdlib>
#include <cstddef>
#include <cstring>
#include <functional>
#include <condition_variable>
#include <chrono>
#include <thread>
#include <cstdint>
#include <algorithm>
#include <cmath>
#include <omp.h>
#include <string>
#include <memory>
#include <limits>
#include <mutex>
#include <numeric>
#include <random>
#include <sstream>
#include <iomanip>


//Random Number Generator
#include "RandomNumberGenerator/lac_param.h"
#include "RandomNumberGenerator/rand.hpp"


//LWE
#include "LWEscheme/LWE_64.hpp"

//RLWE
#include "RLWEscheme/RLWE_64.hpp"


//RLWE
#include "BFVscheme/LT.hpp"
#include "BFVscheme/Poly_eval.hpp"
#include "BFVscheme/LWE_to_RLWE.hpp"
#include "BFVscheme/S2C.hpp"
#include "BFVscheme/S2C_slow.hpp"

//test

#include "test/test_full_protocol.hpp"
#include "test/test_gelu_bert.hpp"
#include "test/test_gelu.hpp"



#endif