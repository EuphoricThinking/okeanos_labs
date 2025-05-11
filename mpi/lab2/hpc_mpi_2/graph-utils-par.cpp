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

static void createGraphFromBuffer(int* buffer, Graph* graph) {//int numVertices, int numProcesses, int srcRank, int numElems) {
    // int firstIdx = getFirstGraphRowOfProcess(numVertices, numProcesses, srcRank);
    // int lastIdxExcl = getFirstGraphRowOfProcess(numVertices, numProcesses, srcRank + 1);

    // Graph* graph = allocateGraphPart(numVertices, firstIdx, lastIdxExcl);

    int numRows = graph->lastRowIdxExcl - graph->firstRowIdxIncl; //numElems / numVertices;
    int numVertices = graph->numVertices;

    // int** rows = (int**) malloc(sizeof(int*) * numRows);
    int idx;
    for (int i = 0; i < numRows; i++) {
        // rows[i] = (int*) malloc(sizeof(int) * numVertices);
        idx = numVertices * i;
        for (int j = 0; j < numVertices; j++) {
            graph->data[i][j] = buffer[idx + j];
        }
    }

    // return graph;
}

static void fillGraphPart(Graph* graph) {
    int numRows = getNumRows(graph);
    int firstIdx = graph->firstRowIdxIncl;
    printf("idx: %d\n", firstIdx);
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
    int numRows = graph->lastRowIdxExcl - graph->firstRowIdxIncl;

    if (myRank == 0) {
        // int numRows = getNumRows(graph);
        // int firstIdx = graph->firstRowIdxIncl;
        // for (int i = 0; i < numRows; ++i) {
        //     initializeGraphRow(graph->data[i], firstIdx, graph->numVertices);
        //     firstIdx++;
        // }
        fillGraphPart(graph);
        int idx = graph->firstRowIdxIncl;
        for (int i = 0; i < numRows; i++) {
            printGraphRow(graph->data[i], idx, numVertices);idx++;
        }
        printf("me\n");
    }

    // int* serializedGraph = getSerializedGraphPart(graph);
    // int numElems = getNumElemsPerGivenProcess(numVertices, numProcesses, myRank);
    int numElems = numRows * numVertices;

    // int maxRowsPerProcess = getMaxNumRowsPerProcess(numVertices, numProcesses);
     //maxRowsPerProcess * numVertices);

    int elemsCountPerGivenProcess;
    int* buffParam;

    Graph* toSend;
    MPI_Status status;
    // the first graph keeps to itself since it is a sender, hence starting from 1
    if (myRank == 0) {
        for (int src = 1; src < numProcesses; ++src) {
            elemsCountPerGivenProcess = getNumElemsPerGivenProcess(numVertices, numProcesses, src);
    
            // if (myRank == 0) {
                toSend = allocateGraphPart(numVertices, getFirstGraphRowOfProcess(numVertices, numProcesses, src), getFirstGraphRowOfProcess(numVertices, numProcesses, src + 1));
                fillGraphPart(toSend);
                buffParam = getSerializedGraphPart(toSend);

                int* iterator = buffParam;
                int srcFirstIdx = getFirstGraphRowOfProcess(numVertices, numProcesses, src);
                for (int i = 0; i < elemsCountPerGivenProcess / numVertices; ++i) {
                    printGraphRow(iterator, srcFirstIdx, numVertices);
                    iterator = iterator + numVertices;
                }
            // }
    
            // buffParam = src == myRank ? serializedGraph : receiveBuffer;
    
            // MPI_Bcast(
            //     receiveBuffer,
            //     elemsCountPerGivenProcess,
            //     MPI_INT,
            //     src,
            //     MPI_COMM_WORLD
            // );
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
        printf("BLAH\n");
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
    // }
    // MPI_Ibcast(
    //     serializedGraph,
    //     numElems,
    //     MPI_INT,
    //     myRank,
    //     MPI_COMM_WORLD
    // // );
    // free(serializedGraph);
    // free(receiveBuffer);



void collectAndPrintGraph(Graph* graph, int numProcesses, int myRank) {
    assert(numProcesses >= 1 && myRank >= 0 && myRank < numProcesses);
    assert(graph->numVertices > 0);
    assert(graph->firstRowIdxIncl >= 0 && graph->lastRowIdxExcl <= graph->numVertices);

    /* FIXME: implement */
    int* serializedGraph = getSerializedGraphPart(graph);
    int numVertices = graph->numVertices;
    int maxRowsPerProcess = getMaxNumRowsPerProcess(numVertices, numProcesses);
    int* receiveBuffer = (int*) malloc(sizeof(int) * maxRowsPerProcess * numVertices);

    int elemsCountPerGivenProcess;
    int* buffParam;
    for (int src = 0; src < numProcesses; ++src) {
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

    free(serializedGraph);
    free(receiveBuffer);
}

void destroyGraph(Graph* graph, int numProcesses, int myRank) {
    freeGraphPart(graph);
}
