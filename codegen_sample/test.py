import dace
from dace import registry
from dace.sdfg.scope import ScopeSubgraphView
from dace.codegen.prettycode import CodeIOStream
from dace.codegen.targets.target import TargetCodeGenerator
from dace.codegen.targets.framecode import DaCeCodeGenerator
from dace.codegen.targets.cpp import sym2cpp
'''
@dace.program
def simple(A: dace.float64[20, 30]):
    for i, j in dace.map[0:20:2, 0:30]:
        A[i, j] += A[i, j]
'''

import numpy as np


m = dace.symbol('m', dace.int32)
P = dace.symbol('P', dace.int32)

@dace.program
def matrix_1d_1d(A: dace.int32[m, m*P, m*P]):
	a_grid = dace.comm.Cart_create([P, 1, 1])
	b_grid = dace.comm.Cart_create([1, P, 1])
	B = np.empty_like(A, shape=(m*P, m, m*P))
	a_arr = dace.comm.Subarray((m*P, m*P, m*P), A, process_grid=a_grid)
	b_arr = dace.comm.Subarray((m*P, m*P, m*P), B, process_grid=b_grid)
	rdistr = dace.comm.Redistribute(A, a_arr, B, b_arr)
	return B
# Preview SDFG
sdfg = matrix_1d_1d.to_sdfg()
from IPython.display import Code
code = Code(sdfg.generate_code()[0].clean_code, language='cpp')
with open('./test3d_3dgrid.cpp', 'w') as f:
    f.write(code.data)
