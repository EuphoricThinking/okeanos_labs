#include "distributed-helpers.h"

int get_num_rows_per_process(int numVertices, int numProcesses) {
    return numVertices / numProcesses;

}

// int min(int a, int b) {
//     if  (a > b) {
//         return b;
//     }

//     return a;
// }

// int max(int a, int b) {
//     if (a > b) {
//         return a;
//     }

//     return b;
// }

int get_last_excl_idx_of_my_part(int numVertices, int numProcesses, int myRank) {
    // the last one always covers the taill
    if (myRank == numProcesses - 1) {
        // exlusive
        return numVertices;
    }
    else {
        int per_process = get_num_rows_per_process(numVertices, numProcesses);

        // int candidate = myRank + per_process;
        // integer division is rounding toward zero, therefore the last process should always cover the remaining tail
        return  myRank + per_process;

        // // our_rank  ver     ver     ver     last_ver   | ...  candidate
        // if (candidate >= numVertices) {
        //     return max(numVertices, myRank + 1);
        // }
        // else {
        //     // if we aer the last process - for corenr caes, when thre is not enough processes to cover all rows
        //     if (numProcesses - 1 == myRank) {
        //         if (candidate < numVertices)  {
        //             return numVertices;
        //         }

        //     }
        // }
    }
}
