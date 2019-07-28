message("\nCMake is only used for building unit tests.\nThe main program is designed for Sunway Cluster and should be built by a specific compiler.")
message("Please contact with chenhaowu@link.cuhk.edu.cn to get more information.\n")

set(CMAKE_CXX_FLAGS "-march=native -Og -fno-stack-protector ")
add_executable(testSolve74 testSolve74.cpp gf2misc.cpp gf2misc.h mqmisc.cpp mqmisc.h mqparam.h)
add_executable(testSolve67 testSolve67.cpp gf2misc.cpp gf2misc.h mqmisc.cpp mqmisc.h mqparam.h)