#ifndef __DISTRIBUTED_HELPERS__H__
#define __DISTRIBUTED_HELPERS__H__

/*
the graph parts are ranndomly enerated between processes



*/

int get_num_rows_per_process(int numVertices, int numProcesses);

int get_last_excl_idx_of_my_part(int numVertices, int numProcesses, int myRank);

#endif /* __DISTRIBUTED_HELPERS__H__ */