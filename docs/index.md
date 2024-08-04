# Overview

NAOMI (<b>N</b>umerical <b>A</b>strodynamics and <b>O</b>rbital <b>M</b>echanics <b>I</b>nterface) is a general purpose library for astrodynamics.
It has been designed from the ground up to support distributed, high fidelity, and interactive
simulations.  Additionally, there is support for symbolic math through the use of [Symengine](https://symengine.org/).

!!! warning

    This project is under active development and subject to significant change without warning.
    It has currently only been tested with gcc on Ubuntu Linux 22.04

## Installation

First install the submodules

```
git submodule init && git submodule update
```

Then install the dependencies.

```
sudo apt-get update
sudo apt-get install libfmt-dev libeigen3-dev libboost-all-dev libopenblas-dev \
    liblapack-dev libarpack2-dev libsuperlu-dev libarmadillo-dev libgtest-dev \
    libgmp-dev libflint-dev libmpc-dev libomp-dev llvm
```

Next build Symengine

```
cd symengine && mkdir build && cd build
cmake WITH_GMP=on WITH_MPFR=on WITH_MPC=on INTEGER_CLASS=flint WITH_LLVM=on WITH_SYMENGINE_THREAD_SAFE=on WITH_OPENMP=on ..
make install
cd ..
```

Finally, build Naomi

```
mkdir build && cd build
cmake ..
make install
cd ..
```

Now the tests can be run from the `build/tests` directory

```
cd build/tests
ctest
```

### Build and Run Examples

To build the examples set `BUILD_EXAMPLES=ON` in the build

```
cmake -DBUILD_EXAMPLES=ON ..
```

In the build directory there will now be an `examples` directory containing
the example executables.  To run the simulation examples:

```
cd build
./examples/simulation/naomi-simulation
```



