#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "topology.h"


#define X_LIM 3
#define Y_LIM 4

#define X1_LIM 2
#define Y1_LIM 3

MPI_Comm* make_top();

int main(int argc, char **argv){
	int i, rank, num_procs;
    MPI_Comm *mesh;

    int err; //FOR ERRORCHECKING

	MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    mesh = make_top();
    if(mesh == NULL){
        //ERROR
        
    }

    int coords[2];

    err = MPI_Cart_coords(*mesh, rank, 2, coords);
    if(err != MPI_SUCCESS){
        printf("MPI_Cart_coords: ERROR\n");
        MPI_Abort(MPI_COMM_WORLD, err);
    }

    err = MPI_Comm_size(*mesh, &num_procs);
    if(err != MPI_SUCCESS){
        printf("MPI_Comm_size: ERROR\n");
        MPI_Abort(MPI_COMM_WORLD, err);
    }


    // printf("%d : (%d,%d)\n", rank, coords[0], coords[1]);

    int temp_rank, cur_coords[2], load_rank[5], coords_ptr[5][2];

    enum CommNum curr_comm;

    MPI_Group mesh_gr, load_gr;
    MPI_Comm load_comm[5];

    /* COMMS FOR ALL */

    for(i=0; i<num_procs; i++){
        MPI_Cart_coords(*mesh, i, 2, cur_coords);
        //printf("%d : (%d,%d) | (%d, %d)\n", rank, coords[0], coords[1], cur_coords[0], cur_coords[1]);
        coords_ptr[LEFT][0]  = cur_coords[0]    ;
        coords_ptr[LEFT][1]  = (cur_coords[1] - 1)%Y_LIM;
        coords_ptr[RIGHT][0] = cur_coords[0]    ;
        coords_ptr[RIGHT][1] = (cur_coords[1] + 1)%Y_LIM;
        coords_ptr[UP][0]    = (cur_coords[0] - 1)%X_LIM;
        coords_ptr[UP][1]    = cur_coords[1]    ;
        coords_ptr[DOWN][0]  = (cur_coords[0] + 1)%X_LIM;
        coords_ptr[DOWN][1]   = cur_coords[1]    ;

        MPI_Cart_rank(*mesh, coords_ptr[LEFT], &load_rank[LEFT]);
        MPI_Cart_rank(*mesh, coords_ptr[RIGHT], &load_rank[RIGHT]);
        MPI_Cart_rank(*mesh, coords_ptr[UP], &load_rank[UP]);
        MPI_Cart_rank(*mesh, coords_ptr[DOWN], &load_rank[DOWN]);

        load_rank[CENTRE] = i;

        if(rank == i){
            //my_color = 1;
            curr_comm = CENTRE;
            //new_rank = CENTRE;
        }
        
        else if(rank == load_rank[LEFT]){
            //my_color = 1;
            curr_comm = RIGHT;
            //new_rank = LEFT;

        }

        else if( rank == load_rank[RIGHT]){
            //my_color = 1;
            curr_comm = LEFT;
            //new_rank = RIGHT;

        }

        else if( rank == load_rank[UP]){
            //my_color = 1;
            curr_comm = DOWN;
            //new_rank = UP;

        }

        else if( rank == load_rank[DOWN]){
            //;
            curr_comm = UP;
            //new_rank = DOWN;

        }

        else {
            //my_color = MPI_UNDEFINED;
            curr_comm = DUMMY;
            //new_rank = DUMMY;
        }

        //printf("%d : (%d,%d) | (%d, %d) [%d]\n", rank, coords[0], coords[1], cur_coords[0], cur_coords[1], my_color);
        //MPI_Comm_split(*mesh, my_color, new_rank, &load_comm[curr_comm]);
        MPI_Comm_group(*mesh, &mesh_gr);
        MPI_Group_incl(mesh_gr, 5, load_rank, &load_gr);
        MPI_Comm_create_group(*mesh, load_gr, 0, &load_comm[curr_comm] );


        if(load_comm[curr_comm] != MPI_COMM_NULL){
            MPI_Comm_rank(load_comm[curr_comm], &temp_rank);
            printf("(%d) %d : %d\n", i, rank, temp_rank);
        }
    }

    int left_rank;
    MPI_Comm_rank(load_comm[UP], &left_rank);
    //printf("%d -> %d\n", rank, left_rank);
    MPI_Finalize();

    return 0;
}


MPI_Comm* make_top(){
    int rank, err;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm *mesh = (MPI_Comm *) malloc(sizeof(MPI_Comm));
    int ndims = 2;
    int dim_sz[2];
    int periodic[2];

    dim_sz[0] = X_LIM;
    dim_sz[1] = Y_LIM;

    periodic[0] = 1;
    periodic[1] = 1;

    int reorder = 1;



    err = MPI_Cart_create(MPI_COMM_WORLD, ndims, dim_sz, periodic, reorder, mesh);
    if(err != MPI_SUCCESS){
        printf("MPI_Cart_create: ERROR\n");
        MPI_Abort(MPI_COMM_WORLD, err);
        return NULL;
        
    }

    return mesh;
}