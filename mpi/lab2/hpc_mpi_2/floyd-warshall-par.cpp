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
#include <stdlib.h>


static int get_relative_idx(int k, Graph* graph) {
    int firstIdx = graph->firstRowIdxIncl;

    return k - firstIdx;
}

static void write_to_k_buffer(int k, Graph* graph, int* k_buffer) {
    int idx = get_relative_idx(k, graph);

    for (int i = 0; i < graph->numVertices; i++) {
        k_buffer[i] = graph->data[idx][i];
    }
}


static void runFloydWarshallParallel(Graph* graph, int numProcesses, int myRank) {
    int numVertices = graph->numVertices;
    assert(numProcesses <= graph->numVertices);


    int numRows = graph->lastRowIdxExcl - graph->firstRowIdxIncl;


    int* kth_row = (int*) malloc(sizeof(int) * numVertices);

    int k_owner;
    int k_owner_row;
    int next_k_owner_row;

    for (k_owner = 0; k_owner < numProcesses; k_owner++) {
        k_owner_row = getFirstGraphRowOfProcess(numVertices, numProcesses, k_owner);
        next_k_owner_row = getFirstGraphRowOfProcess(numVertices, numProcesses, k_owner + 1);

        for (int k = k_owner_row; k < next_k_owner_row; k++) {

            if (k_owner == myRank) {
                write_to_k_buffer(k, graph, kth_row);
            }

            MPI_Bcast(
                kth_row,
                numVertices,
                MPI_INT,
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
