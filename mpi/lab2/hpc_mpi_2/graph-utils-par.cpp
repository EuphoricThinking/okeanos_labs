/*
 * A template for the 2019 MPI lab at the University of Warsaw.
 * Copyright (C) 2016, Konrad Iwanicki.
 * Refactoring 2019, Łukasz Rączkowski
 */

#include <cassert>
#include <mpi.h>
#include "graph-base.h"
#include "graph-utils.h"
#include "distributed-helpers.h"

int getFirstGraphRowOfProcess(int numVertices, int numProcesses, int myRank) {
    /* FIXME: implement */
    auto graph = allocateGraphPart(numVertices, myRank, get_last_excl_idx_of_my_part(numVertices, numProcesses, myRank));

    if (graph == nullptr) {
        return nullptr;
    }

    assert(graph->numVertices > 0 && graph->numVertices == numVertices);

    int my_num_vertices = graph->lastRowIdxExcl - graph->firstRowIdxIncl

    for (int i = graph->firstRowIdxIncl; i < graph->lastRowIdxExcl; ++i) {
        initializeGraphRow(graph->data[i], i, graph->numVertices);
    }


    return myRank;
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

    return graph;
}

void collectAndPrintGraph(Graph* graph, int numProcesses, int myRank) {
    assert(numProcesses >= 1 && myRank >= 0 && myRank < numProcesses);
    assert(graph->numVertices > 0);
    assert(graph->firstRowIdxIncl >= 0 && graph->lastRowIdxExcl <= graph->numVertices);

    /* FIXME: implement */
}

void destroyGraph(Graph* graph, int numProcesses, int myRank) {
    freeGraphPart(graph);
}
