#!/bin/bash

mpicc -o test main.c -lm
./test
rm test
