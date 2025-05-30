
# Overview

This repository provides an interface for applying domain-decomposed solutions of fluid flow ODEs via the Schwarz alternating method through the [pressio-demoapps](https://github.com/Pressio/pressio-demoapps) solver and sample problem suite. This serves as a launching point for exploring Schwarz coupling for advection-dominated systems, as well as coupling full-order ("high-fidelity," FOM) solvers to data-driven projection-based reduced-order models (PROMs) via [Pressio](https://github.com/Pressio/pressio). The framework exemplified in the test cases (in ```tests_cpp/```) should be easily extensible to any sample case provided by **pressio-demoapps**, but as of now is compatible with the 2D shallow water equations, Euler equations, and Burgers' equation. At some point this code may be reworked to apply more generally to codes other than **pressio-demoapps**, but is restricted to this code and cases for now.

# Building and Running Tests

Executing the test cases requires a copy of the **pressio-demoapps** source (which has bundled the **Eigen** library)
and the **pressio**source. Building and executing the test cases can be performed as

```
git clone git@github.com:Pressio/pressio-schwarz.git
export CXX=<path-to-your-CXX-compiler>
export PDA_ROOT=<path-to-pressio-demoapps-root>
export PRESSIO_ROOT=<path-to-pressio>
cd pressio-schwarz && mkdir build && cd build
cmake -DPDA_SOURCE=${PDA_ROOT} -DPRESSIO_SOURCE=${PRESSIO_ROOT} ..
make -j4
ctest -j4
```

# Python utilities

Python utilities for data extraction, visualization, PROM preparation, and error measurement can be found in the ```python/``` directory. Refer to the README there for instructions on installing and using the associated local package.

# Experimental campaign runner

A C++ utility for executing a large number of **pressio-demoapps** and **pressio-schwarz** simulations from YAML input files is provided in the [pdas-experiments](https://github.com/sandialabs/pdas-experiments) repository. This vastly simplifies the numerical experimentation process for the parameterized ODEs supplied by **pressio-demoapps**, and also serves as a repository for input files associated with experimental campaigns presented in publications.