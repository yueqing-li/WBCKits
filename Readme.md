### Requirements
* [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page)
* [qpOASES](https://github.com/coin-or/qpOASES)

### Building and Installation

The WBCKits is build using [CMake](http://www.cmake.org)
To compile the library use:

    mkdir build
    cd build
    cmake ..
    make

To install the library, in the build directory use:

    sudo make install

If can't find the shared library when using, please add the library path like:

    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib

You can test it whether work after installation by using:

    cd test
    mkdir build
    cd build
    cmake ..
    make
    ./xxx_test

### Python Wrapper

If want to use the python wrapper, after install the library use:

    cd python
    python3 setup.py build_ext --inplace

Then copy the pyWBCKits.*.so file to the where you want to use.

You can test it whether work by using:

    python3 test.py
