#include "main.h"
#include "mpi.h"

int main() {
    MPI_Init(0, 0);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    uint64_t segSize = 4503599627370;
    uint64_t init_key = segSize * rank + 2251799813685;
    uint64_t end_key = segSize * (rank + 1);
    mqInit(init_key, end_key);
    mqLoop();
}
