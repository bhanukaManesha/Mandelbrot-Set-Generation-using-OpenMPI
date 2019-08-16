#include <stdio.h>
#include <mpi.h>

int main(int argc, char *argv[]){
    int numtasks, rank;
    int shared_buffer;

    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    MPI_Status status;
    int dim[1], period[1], reorder;
    int coord[2], id;
    MPI_Comm comm;

    int source, dest;

    dim[0]= 4;
    period[0]=1;
    reorder=1;
    MPI_Cart_create(MPI_COMM_WORLD, 1, dim, period, reorder, &comm);

    MPI_Cart_shift( comm, 0, 1, &source, &dest );

	if (source == numtasks - 1) {
        // Is it because of blocking send? the rest hangs for rank 0 to finish?
        scanf( "%d", &shared_buffer );
        MPI_Send(&shared_buffer, 1, MPI_INT, dest, 0, MPI_COMM_WORLD);

    } 
    
    MPI_Recv(&shared_buffer, 1, MPI_INT, source, 0, MPI_COMM_WORLD, &status);
    MPI_Send(&shared_buffer, 1, MPI_INT, dest, 0, MPI_COMM_WORLD);
    
    printf( "Rank %d : %i -> %i -> %i | %d\n", rank, source, rank, dest, shared_buffer);
	// printf( "Rank %d : %d\n", rank, shared_buffer);
    MPI_Finalize();
    return 0;
}

