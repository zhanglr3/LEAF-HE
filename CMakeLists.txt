project(hexl_example LANGUAGES C CXX)
cmake_minimum_required(VERSION 3.13)
set(CMAKE_CXX_STANDARD 14)

include_directories(include)
include_directories(/usr/include)
include_directories(/usr/local/include)
include_directories(/opt/homebrew/Cellar/openssl@3/3.3.0/include)

link_directories(/opt/homebrew/Cellar/openssl@3/3.3.0/lib)
link_directories(/usr/lib)
link_directories(/usr/local/lib) 

set(CMAKE_C_COMPILER "gcc-14")
set(CMAKE_CPP_COMPILER "g++-14")

include(FindPkgConfig)
if(NOT PKG_CONFIG_FOUND)
  message(FATAL_ERROR "pkg-config not found!" )
endif()


#Please add the path of OPENBLAS to the following two lines 
set(CMAKE_CXX_FLAGS "-lstdc++ -O3 -fopenmp ")
set(CMAKE_C_FLAGS "-lstdc++ -O3 -fopenmp ")

#For single-thread test, Please add the path of OPENBLAS to the following two lines
#set(CMAKE_CXX_FLAGS "-lstdc++ -O3  ")
#set(CMAKE_C_FLAGS "-lstdc++ -O3  ")

find_package(SEAL REQUIRED)

add_executable(test test.cpp)
target_link_libraries(test crypto pthread SEAL::seal openblas)



