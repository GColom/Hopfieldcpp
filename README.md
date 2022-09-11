# Hopfield Networks in C++: Requirements and Installation instructions
## Requirements

To use the class implemented in this header it is necessary to have a compiler with support for C++ standard 2011 and OpenMP.

Most compilers come already equipped with these features, in particular no installing of further packages was required during development and testing, on Ubuntu 20.04 using ```g++``` version 9.4.0.

## Installation

To be able to instantiate ```HopfieldNetwork``` objects as implemented in this project, you will need to copy the header file to a directory of your convenience, and ```#include``` it in your programs according to the usual inclusion rules of C++.

Then compile the project with 

```g++ my_program.C -fopenmp -DN_PARALLEL_THREADS=N -Ofast -o my_program.out```

where ```N``` should be the desired default number of parallel threads to use in parallel code sections.
