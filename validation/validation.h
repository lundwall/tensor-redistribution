#pragma once

// This function is only used to validate the receiving value. The sending value could be calculated easily without the help of this header file.
// In the future, we could override this function by a more complicated one once we decided the tensor redistribution task.
int getElement(int i, int rank, int Ni, int Nj, int Ni_new, int Nj_new, int SUBNi, int SUBNj){
    int current_i = i/Nj_new, current_j = i%Nj_new;
    return (current_i < SUBNi && current_j < SUBNj)?current_j+current_i*Nj:0;
}