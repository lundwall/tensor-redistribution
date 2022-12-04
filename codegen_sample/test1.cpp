
#include <mpi.h>
#include <iostream>
#include "dace_helper.h"
        
int main(int argc, char** argv) {
    // we only need change this
    int a0 = 4, a1 = 1, a2 = 1; // pgrid0 dims
    int b0 = 1, b1 = 4, b2 = 1; // pgrid1 dims
    int d0 = 16, d1 = 20, d2 = 24; // data dims

    MPI_Init(&argc, &argv);
    
    MPI_Comm pgrid_0_comm;
    MPI_Comm pgrid_1_comm;

    int pgrid_0_dims[3] = {a0, a1, a2};
    int pgrid_1_dims[3] = {b0, b1, b2};
    int pgrid_0_periods[3] = {0};
    int pgrid_1_periods[3] = {0};
    int pgrid_0_coords[3];
    int pgrid_1_coords[3];
    int pgrid_0_rank;
    int pgrid_1_rank;

    MPI_Cart_create(MPI_COMM_WORLD, 3, pgrid_0_dims, pgrid_0_periods, 0, &pgrid_0_comm);
    MPI_Cart_create(MPI_COMM_WORLD, 3, pgrid_1_dims, pgrid_1_periods, 0, &pgrid_1_comm);
    MPI_Comm_rank(pgrid_0_comm, &pgrid_0_rank);
    MPI_Cart_coords(pgrid_0_comm, pgrid_0_rank, 3, pgrid_0_coords);
    MPI_Comm_rank(pgrid_1_comm, &pgrid_1_rank);
    MPI_Cart_coords(pgrid_1_comm, pgrid_1_rank, 3, pgrid_1_coords);
        
    {
        int kappa[3];
        int lambda[3];
        int xi[3];
        int pcoords[3];
        int myrank;
        MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
        
        {
        
            int sizes[3] = {d0/b0, d1/b1, d2/b2};
            int subsizes[3];
            int origin[3];
        
            xi[0] = (pgrid_1_coords[0] * int(d0/b0)) / int(d0/a0);
            lambda[0] = pgrid_1_coords[0] * int(d0/b0) % int(d0/a0);
            kappa[0] = int_ceil(int(d0/b0) + lambda[0], int(d0/a0));
            
            xi[1] = (pgrid_1_coords[1] * int(d1/b1)) / int(d1/a1);
            lambda[1] = pgrid_1_coords[1] * int(d1/b1) % int(d1/a1);
            kappa[1] = int_ceil(int(d1/b1) + lambda[1], int(d1/a1));
            
            xi[2] = (pgrid_1_coords[2] * int(d2/b2)) / int(d2/a2);
            lambda[2] = pgrid_1_coords[2] * int(d2/b2) % int(d2/a2);
            kappa[2] = int_ceil(int(d2/b2) + lambda[2], int(d2/a2));
            
            int rem0 = d0/b0;
            for (auto idx0 = 0; idx0 < kappa[0]; ++idx0) {
                int actual_idx0 = 0;
                pcoords[actual_idx0] = xi[0] + idx0;
                int lo0 = (idx0 == 0 ? lambda[0] : 0);
                int uo0 = std::min(int(d0/a0), lo0 + rem0);
                subsizes[0] = uo0 - lo0;
                origin[0] = d0/b0 - rem0;
                rem0 -= uo0 - lo0;
            
            int rem1 = d1/b1;
            for (auto idx1 = 0; idx1 < kappa[1]; ++idx1) {
                int actual_idx1 = 1;
                pcoords[actual_idx1] = xi[1] + idx1;
                int lo1 = (idx1 == 0 ? lambda[1] : 0);
                int uo1 = std::min(int(d1/a1), lo1 + rem1);
                subsizes[1] = uo1 - lo1;
                origin[1] = d1/b1 - rem1;
                rem1 -= uo1 - lo1;
            
            int rem2 = d2/b2;
            for (auto idx2 = 0; idx2 < kappa[2]; ++idx2) {
                int actual_idx2 = 2;
                pcoords[actual_idx2] = xi[2] + idx2;
                int lo2 = (idx2 == 0 ? lambda[2] : 0);
                int uo2 = std::min(int(d2/a2), lo2 + rem2);
                subsizes[2] = uo2 - lo2;
                origin[2] = d2/b2 - rem2;
                rem2 -= uo2 - lo2;
            int cart_rank = get_cart_rank(3, pgrid_0_dims, pcoords);

            if (myrank == cart_rank) { // self-copy
                printf("I am rank %d and I self-copy,  receive from %d,%d,%d (%d) from (%d,%d,%d) size (%d,%d,%d)\n", myrank,
                    pcoords[0],pcoords[1],pcoords[2],
                    cart_rank, 
                    origin[0],origin[1],origin[2], 
                    subsizes[0],subsizes[1],subsizes[2]);
            } else {
                printf("I am rank %d and I receive from %d,%d,%d (%d) in (%d,%d,%d) size (%d,%d,%d)\n", myrank,
                    pcoords[0],pcoords[1],pcoords[2],
                    cart_rank, 
                    origin[0],origin[1],origin[2], 
                    subsizes[0],subsizes[1],subsizes[2]);            
            }
        }}}}
        {
        
            int sizes[3] = {d0/a0, d1/a1, d2/a2};
            int subsizes[3];
            int origin[3];
        
            // int_ceil(x, y) := (x + y - 1) / y
            // int_ceil(pcoord * sa - sb + 1, sb) = (pcoord * sa) / sb
            int lp0 = std::max(0, (pgrid_0_coords[0] * int(d0/a0)) / int(d0/b0)); // int_ceil(x, y) := (x + y - 1) / y
            int up0 = std::min(pgrid_1_dims[0], int_ceil((pgrid_0_coords[0] + 1) * int(d0/a0), int(d0/b0)));
            for (auto idx0 = lp0; idx0 < up0; ++idx0) {
                int actual_idx0 = 0;
                xi[0] = (idx0 * int(d0/b0)) / int(d0/a0);
                lambda[0] = idx0 * int(d0/b0) % int(d0/a0);
                kappa[0] = int_ceil(int(d0/b0) + lambda[0], int(d0/a0));
                int idx0_dst = pgrid_0_coords[0] - xi[0];
                if (idx0_dst < 0 || idx0_dst >= kappa[0]) continue;
                int lo0 = (idx0_dst == 0 ? lambda[0] : 0);
                int uo0 = (idx0_dst == kappa[0] - 1 ? int(d0/b0) + lambda[0] - idx0_dst * int(d0/a0) : int(d0/a0));
                subsizes[0] = uo0 - lo0;
                origin[0] = lo0;
                pcoords[actual_idx0] = idx0;
            
            // int_ceil(x, y) := (x + y - 1) / y
            // int_ceil(pcoord * sa - sb + 1, sb) = (pcoord * sa) / sb
            int lp1 = std::max(0, (pgrid_0_coords[1] * int(d1/a1)) / int(d1/b1)); // int_ceil(x, y) := (x + y - 1) / y
            int up1 = std::min(pgrid_1_dims[1], int_ceil((pgrid_0_coords[1] + 1) * int(d1/a1), int(d1/b1)));
            for (auto idx1 = lp1; idx1 < up1; ++idx1) {
                int actual_idx1 = 1;
                xi[1] = (idx1 * int(d1/b1)) / int(d1/a1);
                lambda[1] = idx1 * int(d1/b1) % int(d1/a1);
                kappa[1] = int_ceil(int(d1/b1) + lambda[1], int(d1/a1));
                int idx1_dst = pgrid_0_coords[1] - xi[1];
                if (idx1_dst < 0 || idx1_dst >= kappa[1]) continue;
                int lo1 = (idx1_dst == 0 ? lambda[1] : 0);
                int uo1 = (idx1_dst == kappa[1] - 1 ? int(d1/b1) + lambda[1] - idx1_dst * int(d1/a1) : int(d1/a1));
                subsizes[1] = uo1 - lo1;
                origin[1] = lo1;
                pcoords[actual_idx1] = idx1;
            
            // int_ceil(x, y) := (x + y - 1) / y
            // int_ceil(pcoord * sa - sb + 1, sb) = (pcoord * sa) / sb
            int lp2 = std::max(0, (pgrid_0_coords[2] * int(d2/a2)) / int(d2/b2)); // int_ceil(x, y) := (x + y - 1) / y
            int up2 = std::min(pgrid_1_dims[2], int_ceil((pgrid_0_coords[2] + 1) * int(d2/a2), int(d2/b2)));
            for (auto idx2 = lp2; idx2 < up2; ++idx2) {
                int actual_idx2 = 2;
                xi[2] = (idx2 * int(d2/b2)) / int(d2/a2);
                lambda[2] = idx2 * int(d2/b2) % int(d2/a2);
                kappa[2] = int_ceil(int(d2/b2) + lambda[2], int(d2/a2));
                int idx2_dst = pgrid_0_coords[2] - xi[2];
                if (idx2_dst < 0 || idx2_dst >= kappa[2]) continue;
                int lo2 = (idx2_dst == 0 ? lambda[2] : 0);
                int uo2 = (idx2_dst == kappa[2] - 1 ? int(d2/b2) + lambda[2] - idx2_dst * int(d2/a2) : int(d2/a2));
                subsizes[2] = uo2 - lo2;
                origin[2] = lo2;
                pcoords[actual_idx2] = idx2;
            int cart_rank = get_cart_rank(3, pgrid_1_dims, pcoords);

        if (myrank != cart_rank) { // not self-copy
            printf("I am rank %d and I send to %d,%d,%d (%d) from (%d,%d,%d) size (%d,%d,%d)\n", myrank,
                pcoords[0],pcoords[1],pcoords[2],
                cart_rank, 
                origin[0],origin[1],origin[2], 
                subsizes[0],subsizes[1],subsizes[2]);
        }
        }}}}}   
    MPI_Finalize();
}
        