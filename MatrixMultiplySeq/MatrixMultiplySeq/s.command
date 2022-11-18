#!/bin/bash
mpic++ main.cpp
echo "2 procnum:"
mpiexec -np 2 ./a.out
mpiexec -np 2 ./a.out
mpiexec -np 2 ./a.out
mpiexec -np 2 ./a.out
mpiexec -np 2 ./a.out
echo "4 procnum:"
mpiexec -np 4 ./a.out
mpiexec -np 4 ./a.out
mpiexec -np 4 ./a.out
mpiexec -np 4 ./a.out
mpiexec -np 4 ./a.out
echo "6 procnum:"
mpiexec -np 6 ./a.out
mpiexec -np 6 ./a.out
mpiexec -np 6 ./a.out
mpiexec -np 6 ./a.out
mpiexec -np 6 ./a.out
echo "8 procnum:"
mpiexec -np 8 ./a.out
mpiexec -np 8 ./a.out
mpiexec -np 8 ./a.out
mpiexec -np 8 ./a.out
mpiexec -np 8 ./a.out