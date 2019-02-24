// Copyright 2018 Nesterov Alexander
#include <mpi.h>
#include <stdio.h>

int main(int argc, char* argv[]) {
    int status, rank, size;
    status = MPI_Init(&argc, &argv);
    if (status != MPI_SUCCESS) { return -1; }

    status = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (status != MPI_SUCCESS) { return -1; }

    status = MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (status != MPI_SUCCESS) { return -1; }

    int number;
    if (rank == 0) {
        number = 12;
        MPI_Send(&number, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
    } else if (rank == 1) {
        MPI_Recv(&number, 1, MPI_INT, 0, 0, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
        printf("Process 1 received number %d from process 0\n",
               number);
    }

    status = MPI_Finalize();
    if (status != MPI_SUCCESS) { return -1; }

    return 0;
}
