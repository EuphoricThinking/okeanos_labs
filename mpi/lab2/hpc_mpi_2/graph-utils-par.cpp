/*
 * A template for the 2019 MPI lab at the University of Warsaw.
 * Copyright (C) 2016, Konrad Iwanicki.
 * Refactoring 2019, Łukasz Rączkowski
 */

#include <cassert>
#include <mpi.h>
#include "graph-base.h"
#include "graph-utils.h"
// #include "distributed-helpers.h"


static int getNumRows(Graph* graph) {
    return graph->lastRowIdxExcl - graph->firstRowIdxIncl;
}

static int getMaxNumRowsPerProcess(int numVertices, int numProcesses) {
    int even_per_process = numVertices / numProcesses;
    int rest_vertices = numVertices % numProcesses;

    if (rest_vertices == 0) {
        return even_per_process;
    }
    else {
        return even_per_process + 1;
    }
}

static int* getSerializedGraphPart(Graph* graph) {
    int numRows = getNumRows(graph);
    int numVertices = graph->numVertices;
    int* message = malloc(sizeof(int) * (numVertices * numRows));

    int idx;
    for (int row = 0; row < numRows; ++row) {
        for (int col = 0; col < numVertices; ++col) {
            idx = row * numVertices + col;
            message[idx] = graph->data[row][col];
        }
    }

    // the number of vertices and processes should be known
    return message;
}

static int getNumElemsPerGivenProcess(int numVertices, int numProcesses, int myRank) {
    int even_per_process = numVertices / numProcesses;
    int rest_vertices = numVertices % numProcesses;

    if (rest_vertices == 0 || myRank >= rest_vertices) {
        return even_per_process * numVertices;
    }
    else {
        // if (myRank < rest_vertices) {
            return (even_per_process + 1) * numVertices;
        // }
        // else {
        //     return even_per_process * numVertices;
        // }
    }
}

int getFirstGraphRowOfProcess(int numVertices, int numProcesses, int myRank) {
    /* FIXME: implement */

    // however, main cuts the number of processes to numVertices
    if (numVertices < numProcesses) {
        return myRank < numVertices ? myRank : -1;
    }
    else {
        int even_per_process = numVertices / numProcesses;
        int rest_vertices = numVertices % numProcesses;

        if (rest_vertices == 0) {
            return even_per_process * myRank;
        }
        else {
            // distribute the remaining part evenly among first modulo processes
            if (myRank < rest_vertices) {
                return (even_per_process + 1) * myRank;
            }
            else {
                // the first modulo processes get one more row
                int base = (even_per_process + 1) * rest_vertices;
                int above_modulo = myRank - rest_vertices;

                return base + above_modulo * even_per_process;
            }
        }
    }
    // return per_process * myRank; //myRank;
}

Graph* createAndDistributeGraph(int numVertices, int numProcesses, int myRank) {
    assert(numProcesses >= 1 && myRank >= 0 && myRank < numProcesses);

    auto graph = allocateGraphPart(
            numVertices,
            getFirstGraphRowOfProcess(numVertices, numProcesses, myRank),
            getFirstGraphRowOfProcess(numVertices, numProcesses, myRank + 1)
    );

    if (graph == nullptr) {
        return nullptr;
    }

    assert(graph->numVertices > 0 && graph->numVertices == numVertices);
    assert(graph->firstRowIdxIncl >= 0 && graph->lastRowIdxExcl <= graph->numVertices);

    /* FIXME: implement */
    int numRows = getNumRows(graph);
    int firstIdx = graph->firstRowIdxIncl;
    for (int i = 0; i < numRows; ++i) {
        initializeGraphRow(graph->data[i], firstIdx, graph->numVertices);
        firstIdx++;
    }

    int* serializedGraph = getSerializedGraphPart(graph);
    // int numElems = getNumElemsPerGivenProcess(numVertices, numProcesses, myRank);
    int maxRowsPerProcess = getMaxNumRowsPerProcess(numVertices, numProcesses);
    int* receiveBuffer = malloc(sizeof(int) * maxRowsPerProcess * numVertices);

    int perGivenProcess;
    int* buffParam;
    for (int src = 0; src < numProcesses; ++src {
        elemsCountPerGivenProcess = getNumElemsPerGivenProcess(numVertices, numProcesses, src);
        buffParam = src == myRank ? serializedGraph : receiveBuffer;

        MPI_Bcast(
            buffParam,
            elemsCountPerGivenProcess,
            MPI_INT,
            src,
            MPI_COMM_WORLD
        );
    
    // }
    // MPI_Ibcast(
    //     serializedGraph,
    //     numElems,
    //     MPI_INT,
    //     myRank,
    //     MPI_COMM_WORLD
    // );

    return graph;
}

void collectAndPrintGraph(Graph* graph, int numProcesses, int myRank) {
    assert(numProcesses >= 1 && myRank >= 0 && myRank < numProcesses);
    assert(graph->numVertices > 0);
    assert(graph->firstRowIdxIncl >= 0 && graph->lastRowIdxExcl <= graph->numVertices);

    /* FIXME: implement */
    int* serializedGraph = getSerializedGraphPart(graph);
    int maxRowsPerProcess = getMaxNumRowsPerProcess(numVertices, numProcesses);
    int* receiveBuffer = malloc(sizeof(int) * maxRowsPerProcess * numVertices);

    int perGivenProcess;
    int* buffParam;
    for (int src = 0; src < numProcesses; ++src {
        elemsCountPerGivenProcess = getNumElemsPerGivenProcess(numVertices, numProcesses, src);
        buffParam = src == myRank ? serializedGraph : receiveBuffer;

        MPI_Bcast(
            buffParam,
            elemsCountPerGivenProcess,
            MPI_INT,
            src,
            MPI_COMM_WORLD
        );
    
        if (myRank == 0) {
            int* iterator = buffParam;
            int srcFirstIdx = getFirstGraphRowOfProcess(numVertices, numProcesses, src); // redundant in print Function

            for (int i = 0; i < elemsCountPerGivenProcess / numVertices; ++i) {
                printGraphRow(iterator, srcFirstIdx, numVertices);
                iterator = iterator + numVertices;
            }
        }
    }
}

void destroyGraph(Graph* graph, int numProcesses, int myRank) {
    freeGraphPart(graph);
}
