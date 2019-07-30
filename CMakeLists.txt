message("\nCMake is only used for building unit tests.\nThe main program is designed for Sunway Cluster and should be built by a specific compiler.")
message("Please contact with chenhaowu@link.cuhk.edu.cn to get more information.\n")

set(CMAKE_CXX_FLAGS "-march=native -O -m32 ")
add_executable(testSolve74 SolveType1/testSolve74.cpp SolveType4/gf2misc.cpp SolveType4/gf2misc.h SolveType4/mqmisc.cpp SolveType4/mqmisc.h SolveType4/mqparam.h)
add_executable(testSolve67 SolveType4/testSolve67.cpp SolveType4/gf2misc.cpp SolveType4/gf2misc.h SolveType4/mqmisc.cpp SolveType4/mqmisc.h SolveType4/mqparam.h)
add_executable(solve67Main SolveType4/main.cpp)
add_executable(solve64Main SolveType4_64/main.cpp)