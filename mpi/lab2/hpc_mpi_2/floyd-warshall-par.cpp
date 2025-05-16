/*
 * A template for the 2019 MPI lab at the University of Warsaw.
 * Copyright (C) 2016, Konrad Iwanicki.
 * Refactoring 2019, Łukasz Rączkowski
 */

#include <iostream>
#include <string>
#include <cassert>
#include <mpi.h>
#include <stdbool.h>
#include "graph-base.h"
#include "graph-utils.h"


static bool is_my_row(int row_idx, Graph* graph) {
    return (row_idx >= graph->firstRowIdxIncl) && (row_idx < graph->lastRowIdxExcl);
}

// static void fill_graph_with_inf(Graph* graph) {
//     int num_rows = graph->lastRowIdxExcl - graph->firstRowIdxIncl;

//     for (int row = 0; row < num_rows; row++) {
//         for (int col = 0; col < graph->numVertices; col++) {
//             graph[row][col] = INT_MAX;
//         }
//     }
// }

static int get_relative_idx(int k, Graph* graph) {
    int firstIdx = graph->firstRowIdxIncl;

    return k - firstIdx;
}

static void write_to_k_buffer(int k, Graph* graph, int* k_buffer) {
    int idx = get_relative_idx(k, graph);
    // printf("rel idx: %d num rows %d\n", idx, graph->lastRowIdxExcl - graph->firstRowIdxIncl);

    for (int i = 0; i < graph->numVertices; i++) {
        k_buffer[i] = graph->data[idx][i];
    }
}

static void printer(int* kth_row, int numVertices, int myRank, int k) {
    printf("my rankL %d | k: %d | ", myRank, k);

    for (int i = 0; i < numVertices; i++) {
        printf("%d ", kth_row[i]);
    }

    printf("\n");
}

static int get_owner_of_k(int k, int numVertices, int numProcesses) {
    int per_process = numVertices / numProcesses;
    int rest_processes = numVertices % numProcesses;
    
    if (rest_processes == 0) {
        return k / per_process;
    }
    else {
        int one_more_row_per_process = rest_processes * (per_process + 1);

        if (one_more_row_per_process > k) {
            return k / (per_process + 1);
        }
        else {
            int new_k = k - one_more_row_per_process;
            int base = k / (per_process + 1);
            // int idx = base 
            int res = (new_k / per_process) + base; 
            printf(" k %d >= %d | new_k: %d | base: %d | res: %d\n", k, one_more_row_per_process, new_k, base, res);

            return res;
        }

    }
} 


/*
proc = 4
v = 5

5 / 4 = 1
0 1 | 2 | 3 | 4

k = 4
k_owner = 3
new_k = 4 - 2 = 2;
*/

static void runFloydWarshallParallel(Graph* graph, int numProcesses, int myRank) {
    int numVertices = graph->numVertices;
    assert(numProcesses <= graph->numVertices);

    // auto graph_to_fill = allocateGraphPart(
    //     numVertices,
    //     getFirstGraphRowOfProcess(numVertices, numProcesses, myRank),
    //     getFirstGraphRowOfProcess(numVertices, numProcesses, myRank + 1)
    // );

    // int numRows = graph_to_fill->lastRowIdxExcl - graph_to_fill->firstRowIdxIncl;

    int numRows = graph->lastRowIdxExcl - graph->firstRowIdxIncl;

    // fill_graph_with_inf(graph_to_fill);

    int* kth_row = (int*) malloc(sizeof(int) * numVertices);

    int k_owner;

    for (int k = 0; k < numVertices; k++) {
        k_owner = get_owner_of_k(k, numVertices, numProcesses);
        printf("k owner %d, my rank %d, k %d, my row: %d\n", k_owner, myRank, k,is_my_row(k, graph));


        if (k_owner == myRank) {
        // if (is_my_row(k, graph)) {
            // printf(" is my row: %d rank %d numVer: %d\n", k, myRank, numVertices);
            write_to_k_buffer(k, graph, kth_row);
        }

        MPI_Bcast(
            kth_row,
            numVertices,
            MPI_INT,
            // myRank,
            k_owner,
            MPI_COMM_WORLD
        );

        // printer(kth_row, numVertices, myRank, k);

        for (int i = 0; i < numRows; i++) {
            for (int j = 0; j < numVertices; j++) {
                int pathSum = graph->data[i][k] + kth_row[j];
                if (graph->data[i][j] > pathSum) {
                    graph->data[i][j] = pathSum;
                }
            }
        }
    }
        
    free(kth_row);

    /* FIXME: implement */
}


int main(int argc, char *argv[]) {
    int numVertices = 0;
    int numProcesses = 0;
    int myRank = 0;
    int showResults = 0;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

#ifdef USE_RANDOM_GRAPH
#ifdef USE_RANDOM_SEED
    srand(USE_RANDOM_SEED);
#endif
#endif

    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]).compare("--show-results") == 0) {
            showResults = 1;
        } else {
            numVertices = std::stoi(argv[i]);
        }
    }

    if (numVertices <= 0) {
        std::cerr << "Usage: " << argv[0] << "  [--show-results] <num_vertices>" << std::endl;
        MPI_Finalize();
        return 1;
    }

    if (numProcesses > numVertices) {
        numProcesses = numVertices;

        if (myRank >= numProcesses) {
            MPI_Finalize();
            return 0;
        }
    }

    std::cerr << "Running the Floyd-Warshall algorithm for a graph with " << numVertices << " vertices." << std::endl;

    auto graph = createAndDistributeGraph(numVertices, numProcesses, myRank);
    if (graph == nullptr) {
        std::cerr << "Error distributing the graph for the algorithm." << std::endl;
        MPI_Finalize();
        return 2;
    }

    if (showResults) {
        collectAndPrintGraph(graph, numProcesses, myRank);
    }

    double startTime = MPI_Wtime();

    runFloydWarshallParallel(graph, numProcesses, myRank);

    double endTime = MPI_Wtime();

    std::cerr
            << "The time required for the Floyd-Warshall algorithm on a "
            << numVertices
            << "-node graph with "
            << numProcesses
            << " process(es):\n "
            << endTime - startTime
            << std::endl;

    if (showResults) {
        collectAndPrintGraph(graph, numProcesses, myRank);
    }

    destroyGraph(graph, numProcesses, myRank);

    MPI_Finalize();

    return 0;
}
