# svFSIplus 

<div align="center">

[![Build Status](https://github.com/SimVascular/svFSIplus/actions/workflows/test.yml/badge.svg)](https://github.com/SimVascular/svFSIplus/actions)
[![codecov](https://codecov.io/github/SimVascular/svFSIplus/graph/badge.svg?token=I848DNIHSP)](https://codecov.io/github/SimVascular/svFSIplus)
![Latest Release](https://img.shields.io/github/v/release/SimVascular/svFSIplus?label=latest)
![Platform](https://img.shields.io/badge/platform-macOS%20|%20linux-blue)
[![DOI](https://img.shields.io/badge/DOI-10.21105%2Fjoss.04118-green)](https://doi.org/10.21105/joss.04118)

</div>

## Introduction

svFSIplus is a C++ implementation of the Fortran [svFSI](https://github.com/SimVascular/svFSI) multi-physics finite element solver designed for computational modeling of the cardiovascular system. It represents the <i>First Stage</i> in the development of a C++ multi-physics finite element solver and is essentially a direct line-by-line translation of the [svFSI](https://github.com/SimVascular/svFSI) Fortran code. This was done to maintain a simple mapping between the code of the two versions and to facilitate debugging by comparing intermediate data (i.e. element stiffness matrices) and simulation results.

The *Second Stage* of the solver development will be an entirely new implementation encapsulating portions of the *First Stage* code into C++ classes forming a clean and simple abstraction that can be enhanced by any developer.

## Getting started

Please see the documentation for [getting started](https://simvascular.github.io/svFSIplus/index.html) and [implementation details](https://simvascular.github.io/svFSIplus/implementation.html). To run examples, have a look at our [testing guide](https://simvascular.github.io/svFSIplus/testing.html).

### CMake Options

Three options are available in CMake:

- **buildPy** - Build Python Interface

- **buildDocs** - If ON the documentation is included. After running make, it can be built using
~~~
make docs
~~~

- **sparseSolverType** Use Sparse Matrix Solvers. 

### Sparse Solver Options

The discrete linear system of equations can be solver with various methods. 
These has a significant impact on the solver performance and total solution time. 
The following sparse solver options are available:

- **Skyline Solver**. This is a serial skyline matrix solver. 
- **SuperLU_MT Solver**. This is a sparse solver that performs factorization in parallel on a shared memory machine. ==To use SuperLU_MT you need to download and install its library before building svOneDSolver.==
- **CSparse Solver**. This is a serial sparse solver. The source codes are included with the svOneDSolver source code, so it doesn't require the installation of external libraries. 

#### Skyline Solver

Solutions are typically one order of magniture slower with this solver then other sparse solver options. 

#### SuperLU_MT Solver

SuperLU_MT is a parallel direct linear solver for shared memory architectures. 
The source code for this library is available at  [this page](http://crd-legacy.lbl.gov/~xiaoye/SuperLU/).
The build functionality in svOneDSolver requires the user to specify the local installation folder for the library using the variable **SUPERLU_DIR**. Let's assume you have already downloaded the library, unzipped its content to, say, "path/to/folder/SuperLU_MT_3.0" and built the library through the make utility. From inside the OneDBin folder (see above) you only have to type the following cmake command:

~~~
cmake -DsparseSolverType="superlu" -DSUPERLU_DIR="path/to/folder/SuperLU_MT_3.0" ../svOneDSolver/
~~~

**NOTE**: svOneDSolver assumes that your SuperLU_MT library has been built using the **PTHREAD** library. Please follow the instruction below on how to built SuperLU_MT to make sure this is the case.

==**How to build the SuperLU_MT library so it works with svOneDSolver**==

Please refer to the documentation in the README file distributed with the SuperLU_MT library. 

Here we assume:

- You are in the "path/to/folder/SuperLU_MT_3.0" folder 
- You are using an Ubuntu operating system or equivalent.

For a different operating system and if you want to link SuperLU_MT to a different multithreading library, please refer to the SuperLU_MT installation documentation.

1. **Delete** your file make.inc.
2. **Copy** the file MAKE_INC/make.pthread one folder up and rename it as make.inc.
3. **Install** the BLAS library on your system and edit the BLASLIB entry in make.inc so that it points to the correct location of your system BLAS library:
~~~
BLASLIB = /YOUR-BLAS-PATH/libblas.a
~~~
4. **Build** the code by typing "make blaslib" and then "make". 

#### CSparse

This is the default solver and does not require any external library to be installed. 

CSparse is a coincise sparse linear algebra package. 
The source codes are publicly available at [this link](http://people.sc.fsu.edu/~jburkardt/c_src/csparse/csparse.html).
For convenience, these source code have been included in the svOneDSolver source code.

