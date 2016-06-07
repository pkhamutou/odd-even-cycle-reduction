#!/bin/bash

mpicc -o test parallel.c -lm
mpiexec -n $1 ./test
rm test
