#!/usr/bin/env sh

mpiexec -np 2 build/manual/2dtransmit_manual
mpiexec -np 2 build/manual/3dtransmit_manual
mpiexec -np 2 build/manual/4dtransmit_manual

mpiexec -np 2 build/datatype/2dtransmit
mpiexec -np 2 build/datatype/3dtransmit
mpiexec -np 2 build/datatype/4dtransmit

