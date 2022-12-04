class pgrid(object):
    def __init__(self, name, grid):
        self.name = name
        self.grid = grid


class redistribute_array(object):
    def __init__(self, shape, pgrid, subshape):
        self.shape = shape
        self.pgrid = pgrid
        self.subshape = subshape
        self.correspondence = list(range(len(shape)))


class RedistrArray(object):

    def __init__(self, name, array_a : redistribute_array, array_b : redistribute_array):
        self.name = name
        self.array_a = array_a
        self.array_b = array_b

    def init_code(self):
        array_a = self.array_a
        array_b = self.array_b
        tmp = f"""
#include <mpi.h>
#include <iostream>
#include "dace_helper.h"
        
int main(int argc, char** argv) {{
    // we only need change this
    int a0 = 4, a1 = 1, a2 = 1; // pgrid0 dims
    int b0 = 1, b1 = 4, b2 = 1; // pgrid1 dims
    int d0 = 16, d1 = 20, d2 = 24; // data dims

    MPI_Init(&argc, &argv);
    
    MPI_Comm pgrid_0_comm;
    MPI_Comm pgrid_1_comm;

    int pgrid_0_dims[{len(array_a.pgrid.grid)}] = {{{', '.join([str(s) for s in array_a.pgrid.grid])}}};
    int pgrid_1_dims[{len(array_b.pgrid.grid)}] = {{{', '.join([str(s) for s in array_b.pgrid.grid])}}};
    int pgrid_0_periods[{len(array_a.pgrid.grid)}] = {{0}};
    int pgrid_1_periods[{len(array_b.pgrid.grid)}] = {{0}};
    int pgrid_0_coords[{len(array_a.pgrid.grid)}];
    int pgrid_1_coords[{len(array_b.pgrid.grid)}];
    int pgrid_0_rank;
    int pgrid_1_rank;

    MPI_Cart_create(MPI_COMM_WORLD, {len(array_a.pgrid.grid)}, pgrid_0_dims, pgrid_0_periods, 0, &pgrid_0_comm);
    MPI_Cart_create(MPI_COMM_WORLD, {len(array_b.pgrid.grid)}, pgrid_1_dims, pgrid_1_periods, 0, &pgrid_1_comm);
    MPI_Comm_rank(pgrid_0_comm, &pgrid_0_rank);
    MPI_Cart_coords(pgrid_0_comm, pgrid_0_rank, {len(array_a.pgrid.grid)}, pgrid_0_coords);
    MPI_Comm_rank(pgrid_1_comm, &pgrid_1_rank);
    MPI_Cart_coords(pgrid_1_comm, pgrid_1_rank, {len(array_b.pgrid.grid)}, pgrid_1_coords);
        """
        tmp += f"""
    {{
        int kappa[{len(array_b.shape)}];
        int lambda[{len(array_b.shape)}];
        int xi[{len(array_b.shape)}];
        int pcoords[{len(array_b.shape)}];
        int myrank;
        MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
        """
        tmp += f"""
        {{
        """
        tmp += f"""
            int sizes[{len(array_b.subshape)}] = {{{', '.join([str(s) for s in array_b.subshape])}}};
            int subsizes[{len(array_b.subshape)}];
            int origin[{len(array_b.subshape)}];
        """
        for i, (sa, sb, cb) in enumerate(zip(array_a.subshape, array_b.subshape, array_b.correspondence)):
            sa = f"int({sa})"
            sb = f"int({sb})"
            pcoord = f"{array_b.pgrid.name}_coords[{cb}]"
            tmp += f"""
            xi[{i}] = ({pcoord} * {sb}) / {sa};
            lambda[{i}] = {pcoord} * {sb} % {sa};
            kappa[{i}] = int_ceil({sb} + lambda[{i}], {sa});
            """
        for i in range(len(array_b.shape)):
            tmp += f"""
            int rem{i} = {array_b.subshape[i]};
            for (auto idx{i} = 0; idx{i} < kappa[{i}]; ++idx{i}) {{
                int actual_idx{i} = {array_a.correspondence[i]};
                pcoords[actual_idx{i}] = xi[{i}] + idx{i};
                int lo{i} = (idx{i} == 0 ? lambda[{i}] : 0);
                int uo{i} = std::min(int({array_a.subshape[i]}), lo{i} + rem{i});
                subsizes[{i}] = uo{i} - lo{i};
                origin[{i}] = {array_b.subshape[i]} - rem{i};
                rem{i} -= uo{i} - lo{i};
            """

        grid_a = array_a.pgrid.grid
        tmp += f"int cart_rank = get_cart_rank({len(grid_a)}, {array_a.pgrid.name}_dims, pcoords);\n"
        tmp += f"""
            if (myrank == cart_rank) {{ // self-copy
                printf("I am rank %d and I self-copy,  receive from {','.join(['%d'for i in range(len(array_a.pgrid.grid))])} (%d) from ({','.join(['%d'for i in range(len(array_a.shape))])}) size ({','.join(['%d'for i in range(len(array_a.shape))])})\\n", myrank,
                    {','.join([('pcoords['+str(i)+']') for i in range(len(array_a.pgrid.grid))])},
                    cart_rank, 
                    {','.join([('origin['+str(i)+']') for i in range(len(array_a.shape))])}, 
                    {','.join([('subsizes['+str(i)+']') for i in range(len(array_a.shape))])});
            }} else {{
                printf("I am rank %d and I receive from {','.join(['%d'for i in range(len(array_a.pgrid.grid))])} (%d) in ({','.join(['%d'for i in range(len(array_a.shape))])}) size ({','.join(['%d'for i in range(len(array_a.shape))])})\\n", myrank,
                    {','.join([('pcoords['+str(i)+']') for i in range(len(array_a.pgrid.grid))])},
                    cart_rank, 
                    {','.join([('origin['+str(i)+']') for i in range(len(array_a.shape))])}, 
                    {','.join([('subsizes['+str(i)+']') for i in range(len(array_a.shape))])});            
            }}
        """
        for i in range(len(array_b.shape)):
            tmp += f"}}"
        tmp += "}"
        tmp += f"""
        {{
        """
        tmp += f"""
            int sizes[{len(array_a.subshape)}] = {{{', '.join([str(s) for s in array_a.subshape])}}};
            int subsizes[{len(array_a.subshape)}];
            int origin[{len(array_a.subshape)}];
        """
        for i in range(len(array_b.shape)):
            pcoord = f"{array_a.pgrid.name}_coords[{array_a.correspondence[i]}]"
            sa = f"int({array_a.subshape[i]})"
            sb = f"int({array_b.subshape[i]})"
            tmp += f"""
            // int_ceil(x, y) := (x + y - 1) / y
            // int_ceil(pcoord * sa - sb + 1, sb) = (pcoord * sa) / sb
            int lp{i} = std::max(0, ({pcoord} * {sa}) / {sb}); // int_ceil(x, y) := (x + y - 1) / y
            int up{i} = std::min({array_b.pgrid.name}_dims[{array_b.correspondence[i]}], int_ceil(({pcoord} + 1) * {sa}, {sb}));
            for (auto idx{i} = lp{i}; idx{i} < up{i}; ++idx{i}) {{
                int actual_idx{i} = {array_b.correspondence[i]};
                xi[{i}] = (idx{i} * {sb}) / {sa};
                lambda[{i}] = idx{i} * {sb} % {sa};
                kappa[{i}] = int_ceil({sb} + lambda[{i}], {sa});
                int idx{i}_dst = {pcoord} - xi[{i}];
                if (idx{i}_dst < 0 || idx{i}_dst >= kappa[{i}]) continue;
                int lo{i} = (idx{i}_dst == 0 ? lambda[{i}] : 0);
                int uo{i} = (idx{i}_dst == kappa[{i}] - 1 ? {sb} + lambda[{i}] - idx{i}_dst * {sa} : {sa});
                subsizes[{i}] = uo{i} - lo{i};
                origin[{i}] = lo{i};
                pcoords[actual_idx{i}] = idx{i};
            """

        grid_b = array_b.pgrid.grid
        tmp += f"int cart_rank = get_cart_rank({len(grid_b)}, {array_b.pgrid.name}_dims, pcoords);\n"
        tmp += f"""
        if (myrank != cart_rank) {{ // not self-copy
            printf("I am rank %d and I send to {','.join(['%d'for i in range(len(array_b.pgrid.grid))])} (%d) from ({','.join(['%d'for i in range(len(array_b.shape))])}) size ({','.join(['%d'for i in range(len(array_b.shape))])})\\n", myrank,
                {','.join([('pcoords['+str(i)+']') for i in range(len(array_b.pgrid.grid))])},
                cart_rank, 
                {','.join([('origin['+str(i)+']') for i in range(len(array_b.shape))])}, 
                {','.join([('subsizes['+str(i)+']') for i in range(len(array_b.shape))])});
        }}
        """
        for i in range(len(array_b.shape)):
            tmp += f"}}"
        tmp += "}"
        tmp += "}"
        tmp += f"""   
    MPI_Finalize();
}}
        """
        return tmp


name = "__rdistrarray_0"


num_dim = 3 # we make the dim of data and pgrid equal
# array_a = redistribute_array(["d0", "d1", "d2"], pgrid("pgrid_0", ["a0", "a1", "a2"]), ["d0/a0", "d1/a1", "d2/a2"])
# array_b = redistribute_array(["d0", "d1", "d2"], pgrid("pgrid_1", ["b0", "b1", "b2"]), ["d0/b0", "d1/b1", "d2/b2"])
array_a = redistribute_array([('d'+str(i)) for i in range(num_dim)], pgrid("pgrid_0", [('a'+str(i)) for i in range(num_dim)]), [('d'+str(i)+'/a'+str(i)) for i in range(num_dim)])
array_b = redistribute_array([('d'+str(i)) for i in range(num_dim)], pgrid("pgrid_1", [('b'+str(i)) for i in range(num_dim)]), [('d'+str(i)+'/b'+str(i)) for i in range(num_dim)])


t = RedistrArray(name, array_a, array_b)
print(t.init_code())
with open('./test1.cpp', 'w') as f:
    f.write(t.init_code())