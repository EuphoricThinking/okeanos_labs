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
            // distribute the remaining part evenly among modulo first processes
            if (myRank < rest_vertices) {
                return (even_per_process + 1) * myRank;
            }
            else {
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
    for (int i = 0; i < graph->numVertices; ++i) {
        initializeGraphRow(graph->data[i], i, graph->numVertices);
    }

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
