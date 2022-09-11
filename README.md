# Hopfield Networks in C++: Requirements and Installation instructions
## Requirements

To use the class implemented in this header it is necessary to have a compiler with support for C++ standard 2011 and OpenMP.

Most compilers come already equipped with these features, in particular no installing of further packages was required during development and testing, on Ubuntu 20.04 using ```g++``` version 9.4.0.

## Installation

To be able to instantiate ```HopfieldNetwork``` objects as implemented in this project, you will need to copy the header file to a directory of your convenience, and ```#include``` it in your programs according to the usual inclusion rules of C++.

Then compile the project with 

```g++ my_program.C -fopenmp -DN_PARALLEL_THREADS=N -Ofast -o my_program.out```

where ```N``` should be the desired default number of parallel threads to use in parallel code sections.

## Running the example and tests

To run the example and the tests, please create a `Hopfieldcpp/build` directory first.
The example program and the tests can be compiled by navigating to `Hopfieldcpp/source/` and launching 

`cmake -D N_PARALLEL_THREADS=N CMakeLists.txt`

where ```N``` should be the desired default number of parallel threads to use in parallel code sections (`N` defaults to the result of the `nproc` command in absence of the `-D` option). 

The output of `cmake` should follow, after that pass `make` to build the example or `make test` to run the automated tests.

The example program can be run by passing `./build/example` from the `Hopfieldcpp` directory.

The example program consists of a program that calculates the pure state amplitude graph, it has the following calling signature

`./example N_spins alpha T_lower_bound T_upper_bound T_steps repetitions_per_step`

and produces a table to standard output.
