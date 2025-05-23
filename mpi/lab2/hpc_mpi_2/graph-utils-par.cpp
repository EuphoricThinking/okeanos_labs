/*
 * A template for the 2019 MPI lab at the University of Warsaw.
 * Copyright (C) 2016, Konrad Iwanicki.
 * Refactoring 2019, Łukasz Rączkowski
 */

#include <cassert>
#include <mpi.h>
#include "graph-base.h"
#include "graph-utils.h"
#include <stdlib.h>


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
    int* message = (int*) malloc(sizeof(int) * (numVertices * numRows));

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

static void createGraphFromBuffer(int* buffer, Graph* graph) {
    int numRows = graph->lastRowIdxExcl - graph->firstRowIdxIncl;
    int numVertices = graph->numVertices;

    int idx;
    for (int i = 0; i < numRows; i++) {
        idx = numVertices * i;
        for (int j = 0; j < numVertices; j++) {
            graph->data[i][j] = buffer[idx + j];
        }
    }
}

static void fillGraphPart(Graph* graph) {
    int numRows = getNumRows(graph);
    int firstIdx = graph->firstRowIdxIncl;
    for (int i = 0; i < numRows; ++i) {
        initializeGraphRow(graph->data[i], firstIdx, graph->numVertices);
        firstIdx++;
    }
}

static int getNumElemsPerGivenProcess(int numVertices, int numProcesses, int myRank) {
    int even_per_process = numVertices / numProcesses;
    int rest_vertices = numVertices % numProcesses;

    if (rest_vertices == 0 || myRank >= rest_vertices) {
        return even_per_process * numVertices;
    }
    else {
        return (even_per_process + 1) * numVertices;
    }
}

static int min(int a, int b) {
    if (a > b) return b;
    return a;
}

int getFirstGraphRowOfProcess(int numVertices, int numProcesses, int myRank) {
    /* FIXME: implement */
    int even_per_process = numVertices / numProcesses;
    int rest_vertices = numVertices % numProcesses;

    return myRank * even_per_process + min(myRank, rest_vertices);
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
    int numRows = graph->lastRowIdxExcl - graph->firstRowIdxIncl;

    if (myRank == 0) {
        fillGraphPart(graph);
    }

    int numElems = numRows * numVertices;

    int elemsCountPerGivenProcess;
    int* buffParam;

    Graph* toSend;
    MPI_Status status;
    // the first graph keeps to itself since it is a sender, hence starting from 1
    if (myRank == 0) {
        for (int src = 1; src < numProcesses; ++src) {
            elemsCountPerGivenProcess = getNumElemsPerGivenProcess(numVertices, numProcesses, src);
    
            toSend = allocateGraphPart(numVertices, getFirstGraphRowOfProcess(numVertices, numProcesses, src), getFirstGraphRowOfProcess(numVertices, numProcesses, src + 1));
            fillGraphPart(toSend);
            buffParam = getSerializedGraphPart(toSend);

            MPI_Send(buffParam,  /* pointer to the message */
                elemsCountPerGivenProcess, /* number of items in the message */
                MPI_INT, /* type of data in the message */
                src, /* rank of the destination process */
                0, /* app-defined message type */
                MPI_COMM_WORLD /* communicator to use */
            );

            freeGraphPart(toSend);
            free(buffParam);
        }
    }
    else {
        int* receiveBuffer = (int*) malloc(sizeof(int) * numElems);
         MPI_Recv(receiveBuffer, /* where the message will be saved */
            numElems, /* max number of elements we expect */
            MPI_INT, /* type of data in the message */
            0, /* if not MPI_ANY_SOURCE, receive only from source with the given rank  */
            0, /* if not MPI_ANY_TAG, receive only with a certain tag */
            MPI_COMM_WORLD, /* communicator to use */
            &status /* if not MPI_STATUS_IGNORE, write comm info here */
        );

        createGraphFromBuffer(receiveBuffer, graph);
        free(receiveBuffer);
    }

    return graph;
}

void collectAndPrintGraph(Graph* graph, int numProcesses, int myRank) {
    assert(numProcesses >= 1 && myRank >= 0 && myRank < numProcesses);
    assert(graph->numVertices > 0);
    assert(graph->firstRowIdxIncl >= 0 && graph->lastRowIdxExcl <= graph->numVertices);

    /* FIXME: implement */
    int numVertices = graph->numVertices;
    MPI_Status status;
    int numRows = graph->lastRowIdxExcl - graph->firstRowIdxIncl;

    if (myRank == 0) {
        int maxRowsPerProcess = getMaxNumRowsPerProcess(numVertices, numProcesses);
        int* receiveBuffer = (int*) malloc(sizeof(int) * maxRowsPerProcess * numVertices);
        int elemsCountPerGivenProcess;

        // first: print myself
        int myIdx = graph->firstRowIdxIncl;
        for (int i = 0; i < numRows; ++i) {
            printGraphRow(graph->data[i], myIdx, graph->numVertices);
            myIdx++;
        }

        for (int src = 1; src < numProcesses; ++src) {
            elemsCountPerGivenProcess = getNumElemsPerGivenProcess(numVertices, numProcesses, src);
         MPI_Recv(receiveBuffer, /* where the message will be saved */
            elemsCountPerGivenProcess, /* max number of elements we expect */
            MPI_INT, /* type of data in the message */
            src, /* if not MPI_ANY_SOURCE, receive only from source with the given rank  */
            0, /* if not MPI_ANY_TAG, receive only with a certain tag */
            MPI_COMM_WORLD, /* communicator to use */
            &status /* if not MPI_STATUS_IGNORE, write comm info here */

           );

            int* iterator = receiveBuffer; //buffParam;
            int srcFirstIdx = getFirstGraphRowOfProcess(numVertices, numProcesses, src); // redundant in print Function

            for (int i = 0; i < elemsCountPerGivenProcess / numVertices; ++i) {
                printGraphRow(iterator, srcFirstIdx, numVertices);
                iterator = iterator + numVertices;
            }
        }

        free(receiveBuffer);
    }
    else {
        int* serializedGraph = getSerializedGraphPart(graph);
        int numElems = numRows * numVertices;

        MPI_Send(serializedGraph,  /* pointer to the message */
            numElems, /* number of items in the message */
            MPI_INT, /* type of data in the message */
            0, /* rank of the destination process */
            0, /* app-defined message type */
            MPI_COMM_WORLD /* communicator to use */
           );

           free(serializedGraph);
    }
}

void destroyGraph(Graph* graph, int numProcesses, int myRank) {
    freeGraphPart(graph);
}
