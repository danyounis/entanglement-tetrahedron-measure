#!/bin/bash
# depends: anaconda3/2021.11, gcc/11.2.0/b1, lapack/3.9.0/b2, openblas/0.3.10/b1
python -m numpy.f2py -llapack -lblas -c Primitives.f95 -m Primitives
