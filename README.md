# NAOMI

sudo apt-get install libfmt-dev libeigen3-dev libboost-all-dev libopenblas-dev 
    liblapack-dev libarpack2-dev libsuperlu-dev libarmadillo-dev

### Symengine

sudo apt-get install libgmp-dev libflint-dev libmpc-dev libomp-dev llvm

    mkdir build && cd build
    cmake -DWITH_GMP=on -DWITH_MPFR=on -DWITH_MPC=on -DINTEGER_CLASS=flint -DWITH_LLVM=on -DWITH_SYMENGINE_THREAD_SAFE=on -DWITH_OPENMP=on ..    make
    make install