#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include <omp.h>
#include <mpi.h>

#include "priority_queue.h"
#include "generics.h"
#include "comforts.h"

#define SYN_TAG 11
#define SYN_ACK_TAG 22
#define DATA_TAG 33
#define INIT_LOAD_TAG 44

#define NUM_THREADS 9

int working[NUM_THREADS];
int can_work[NUM_THREADS];
queue* shared_queue;

int err = 0;
int my_rank;
int world_size;
int X_LIM;
int Y_LIM;
int end_program = 0;

void queue_balancing(int neighbor_rank, int my_queue_len, int neighbour_queue_len);
void expand_partial_solution(queue* private_queue, void* domain_data);
MPI_Comm* create_topology();
void setupCommunicators();
void load_balancing();
int initializeLoad(void* domain_data);

enum CommNum {
	CENTRE,
	LEFT,
	RIGHT,
	UP,
	DOWN,
	DUMMY
};

typedef struct bound_comm {
	float new_bound;
	int start_rank;
} bound_comm;

MPI_Datatype bound_comm_t;

omp_lock_t best_solution_lock;
int best_score_changed = 0;
float best_score = FLT_MAX;
solution_vector best_solution;
MPI_Comm torus;
int neighbors_rank[5];
int stop_flag = 0;


int main(int argc, char** argv) {
	MPI_Init(&argc, &argv);

	setupCommunicators();
	MPI_Comm_rank(torus, &my_rank);
	srand(time(NULL));

	omp_init_lock(&best_solution_lock);

	void* domain_data = populate_domain_data(argc, argv);

	shared_queue = create_queue();
	for (int i = 0; i < NUM_THREADS; i++) {
		working[i] = 0;
		can_work[i] = 0;
	}

	initializeLoad(domain_data);
// 	load_balancing();

	MPI_Finalize();

	return 0;
}

void setupCommunicators() {

	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	X_LIM = 1;
	Y_LIM = world_size;

	for (int i = ((int) sqrt(world_size)); i > 1; i--) {
		if (world_size % i == 0) {
			X_LIM = i;
			Y_LIM = world_size/X_LIM;
			break;
		}
	}

	int rank;
	torus = *(create_topology());
	MPI_Comm_rank(torus, &rank);
	int coords[2], cur_coords[2];

	REPORT_ERROR(MPI_Cart_coords(torus, rank, 2, coords), "MPI_Cart_coords: ERROR\n", );
	int coords_ptr[5][2];

	coords_ptr[CENTRE][0] = coords[0];
	coords_ptr[CENTRE][1] = coords[1];
	coords_ptr[LEFT][0]  = coords[0];
	coords_ptr[LEFT][1]  = (coords[1] - 1);
	coords_ptr[RIGHT][0] = coords[0];
	coords_ptr[RIGHT][1] = (coords[1] + 1);
	coords_ptr[UP][0]    = (coords[0] - 1);
	coords_ptr[UP][1]    = coords[1];
	coords_ptr[DOWN][0]  = (coords[0] + 1);
	coords_ptr[DOWN][1]   = coords[1];

	for (int i = 0; i < 5; i++) {
		if ((coords_ptr[i][0] >= 0 && coords_ptr[i][0] < X_LIM) &&
			(coords_ptr[i][1] >= 0 && coords_ptr[i][1] < Y_LIM)) {
			MPI_Cart_rank(torus, coords_ptr[i], &neighbors_rank[i]);
		} else {
			neighbors_rank[i] = -1;
		}
	}
}

MPI_Comm* create_topology() {
	MPI_Comm *mesh = (MPI_Comm *) malloc(sizeof(MPI_Comm));
	int ndims = 2;
	int dim_sz[2];
	int periodic[2];

	dim_sz[0] = X_LIM;
	dim_sz[1] = Y_LIM;

	periodic[0] = 0;
	periodic[1] = 0;

	int reorder = 1;

	REPORT_ERROR(MPI_Cart_create(MPI_COMM_WORLD, ndims, dim_sz, periodic, reorder, mesh),
				 "MPI_Cart_create: ERROR\n",
			  NULL);

	return mesh;
}

void expand_partial_solution(queue*  private_queue, void* domain_data) {
	float score;
	solution_vector partial_solution = pq_min_extract(shared_queue, &score);

	if (partial_solution != NULL) {
		if (construct_candidates(partial_solution, score, private_queue, domain_data)) {
			int flag = 0;

		omp_set_lock(&best_solution_lock);
			if (score < best_score) {
				best_score = score;
				best_solution = partial_solution;
				best_score_changed++;
			flag = 1;
		}
		omp_unset_lock(&best_solution_lock);

			if (flag) {
				pq_prune(shared_queue, score);
			}
		} else {
			pq_prune(private_queue, best_score);
			pq_merge(shared_queue, private_queue);
		}
	}
}

int getQueueLen() {
	return rand()%100;
}

// void random_perm(int* order, int len) {
// 	int temp;
// 	int a, b;
// 	for (int i = 0; i < len; i++) {
// 		a = rand()%len;
// 		b = rand()%len;
// 		SWAP(order[a], order[b]);
// 	}
// }

// void load_balancing() {
// 	int syn_share = 1;
// 	int syn_receive;
// 	int syn_ack_receive;
// 	int syn_ack_send = 1;
// 	int recv_queue_len;
// 	int send_queue_len;
// 	int flag = 0;
//
// 	MPI_Request syn_request_right, syn_request_down;
// 	MPI_Request syn_ack_right, syn_ack_down;
// 	MPI_Request syn_receive_left, syn_receive_up;
//
// 	if (neighbors_rank[RIGHT] != -1) {
// 		MPI_Isend(&syn_share, 1, MPI_INT, neighbors_rank[RIGHT], SYN_TAG, torus,
// 					&syn_request_right);
// 	}
// 	if (neighbors_rank[DOWN] != -1) {
// 		MPI_Isend(&syn_share, 1, MPI_INT, neighbors_rank[DOWN], SYN_TAG,
// 				torus, &syn_request_down);
// 	}
//
// 	if (neighbors_rank[LEFT] != -1) {
// 		MPI_Irecv(&syn_receive, 1, MPI_INT, neighbors_rank[LEFT], SYN_TAG, torus,
// 					&syn_receive_left);
// 	}
//
// 	if (neighbors_rank[UP] != -1) {
// 		MPI_Irecv(&syn_receive, 1, MPI_INT, neighbors_rank[UP], SYN_TAG, torus,
// 				  &syn_receive_up);
// 	}
//
// 	if (neighbors_rank[RIGHT] != -1) {
// 		MPI_Irecv(&syn_ack_receive, 1, MPI_INT, neighbors_rank[RIGHT], SYN_ACK_TAG, torus,
// 				  &syn_ack_right);
// 	}
//
// 	if (neighbors_rank[DOWN] != -1) {
// 		MPI_Irecv(&syn_ack_receive, 1, MPI_INT, neighbors_rank[DOWN], SYN_ACK_TAG, torus,
// 				  &syn_ack_down);
// 	}
//
// 	while (!stop_flag) {
// 		if (neighbors_rank[LEFT] != -1) {
// 			MPI_Test(&syn_receive_left, &flag, MPI_STATUS_IGNORE);
//
// 			if (flag == 1) {
// 				MPI_Send(&syn_ack_send, 1, MPI_INT, neighbors_rank[LEFT], SYN_ACK_TAG, torus);
// 				MPI_Recv(&recv_queue_len, 1, MPI_INT, neighbors_rank[LEFT], DATA_TAG, torus,
// 						 MPI_STATUS_IGNORE);
// 				send_queue_len = getQueueLen();
// 				MPI_Send(&send_queue_len, 1, MPI_INT, neighbors_rank[LEFT], DATA_TAG, torus);
// 				queue_balancing(neighbors_rank[LEFT], send_queue_len, recv_queue_len);
//
// 				MPI_Irecv(&syn_receive, 1, MPI_INT, neighbors_rank[LEFT], SYN_TAG,
// 						  torus, &syn_receive_left);
// 			}
// 		}
//
// 		if (neighbors_rank[UP] != -1) {
// 			MPI_Test(&syn_receive_up, &flag, MPI_STATUS_IGNORE);
//
// 			if (flag) {
// 				MPI_Send(&syn_ack_send, 1, MPI_INT, neighbors_rank[UP], SYN_ACK_TAG, torus);
// 				MPI_Recv(&recv_queue_len, 1, MPI_INT, neighbors_rank[UP], DATA_TAG, torus,
// 						 MPI_STATUS_IGNORE);
// 				send_queue_len = getQueueLen();
// 				MPI_Send(&send_queue_len, 1, MPI_INT, neighbors_rank[UP], DATA_TAG, torus);
// 				queue_balancing(neighbors_rank[UP], send_queue_len, recv_queue_len);
//
// 				MPI_Irecv(&syn_receive, 1, MPI_INT, neighbors_rank[UP], SYN_TAG,
// 						  torus, &syn_receive_up);
// 			}
// 		}
//
// 		if (neighbors_rank[RIGHT] != -1) {
// 			MPI_Test(&syn_ack_right, &flag, MPI_STATUS_IGNORE);
//
// 			if (flag) {
// 				send_queue_len = getQueueLen();
// 				MPI_Send(&send_queue_len, 1, MPI_INT, neighbors_rank[RIGHT], DATA_TAG, torus);
// 				MPI_Recv(&recv_queue_len, 1, MPI_INT, neighbors_rank[RIGHT], DATA_TAG, torus,
// 						 MPI_STATUS_IGNORE);
//
// 				queue_balancing(neighbors_rank[RIGHT], send_queue_len, recv_queue_len);
//
// 				MPI_Isend(&syn_share, 1, MPI_INT, neighbors_rank[RIGHT], SYN_TAG, torus,
// 						&syn_request_right);
// 				MPI_Irecv(&syn_ack_receive, 1, MPI_INT, neighbors_rank[RIGHT], SYN_ACK_TAG,
// 						  torus, &syn_ack_right);
// 			}
// 		}
//
// 		if (neighbors_rank[DOWN] != -1) {
// 			MPI_Test(&syn_ack_down, &flag, MPI_STATUS_IGNORE);
//
// 			if (flag) {
// 				send_queue_len = getQueueLen();
// 				MPI_Send(&send_queue_len, 1, MPI_INT, neighbors_rank[DOWN], DATA_TAG, torus);
// 				MPI_Recv(&recv_queue_len, 1, MPI_INT, neighbors_rank[DOWN], DATA_TAG, torus,
// 						 MPI_STATUS_IGNORE);
//
// 				queue_balancing(neighbors_rank[DOWN], send_queue_len, recv_queue_len);
//
// 				MPI_Isend(&syn_share, 1, MPI_INT, neighbors_rank[DOWN], SYN_TAG,
// 						  torus, &syn_request_down);
// 				MPI_Irecv(&syn_ack_receive, 1, MPI_INT, neighbors_rank[DOWN], SYN_ACK_TAG,
// 						  torus, &syn_ack_down);
// 			}
// 		}
// 	}
// }
//
//
// void queue_balancing(int neighbor_rank, int my_queue_len, int neighbour_queue_len) {
// 	printf("RANK:%d NEIGHBOUR_RANK:%d MY_QUEUE:%d NEIGHBOUR_QUEUE:%d \n", my_rank,
// 			neighbor_rank, my_queue_len, neighbour_queue_len);
//
// }

// void send_bound() {
// 	MPI_Request request;
// 	bound_comm recv_bound;
// 	recv_bound->start_rank = my_rank;
//
// 	omp_set_lock(&best_solution_lock);
// 	recv_bound->new_bound = best_score;
// 	best_score_changed = 0;
// 	omp_unset_lock(&best_solution_lock);
//
// 	REPORT_ERROR(MPI_Isend(&recv_bound, 1, bound_comm_t,  getUpRank(my_rank, torus),
// 				BOUND_TAG, torus, &request), "MPI_Isend : ERROR\n", );
// 	REPORT_ERROR(MPI_Isend(&recv_bound, 1, bound_comm_t,  getDownRank(my_rank, torus),
// 						   BOUND_TAG, torus, &request), "MPI_Isend : ERROR\n", );
// 	REPORT_ERROR(MPI_Isend(&recv_bound, 1, bound_comm_t,  getLeftRank(my_rank, torus),
// 						   BOUND_TAG, torus, &request), "MPI_Isend : ERROR\n", );
// 	REPORT_ERROR(MPI_Isend(&recv_bound, 1, bound_comm_t,  getRightRank(my_rank, torus),
// 						   BOUND_TAG, torus, &request), "MPI_Isend : ERROR\n", );
// }
//
// void receive_and_forward_bound() {
// 	int recv_flag = 0;
// 	int flag = 0;
// 	MPI_Status status;
// 	MPI_Request request;
// 	bound_comm recv_bound;
//
// 	REPORT_ERROR(MPI_Iprobe(getLeftRank(my_rank, torus), BOUND_TAG, torus, &recv_flag,
// 							&status), "MPI_Iprobe : ERROR\n", );
// 	if (recv_flag) {
// 		REPORT_ERROR(MPI_Recv(&recv_bound, 1, bound_comm_t, getLeftRank(my_rank, torus),
// 							BOUND_TAG, torus, MPI_STATUS_IGNORE),
// 					"MPI_Recv : ERROR\n", );
//
// 		flag = 0;
// 		omp_set_lock(&best_solution_lock);
// 		if (recv_bound->new_bound < best_score) {
// 			best_score = recv_bound->new_bound;
// 			best_solution = NULL;
// 			flag = 1;
// 		}
// 		omp_unset_lock(&best_solution_lock);
//
// 		if (flag == 1) {
// 			pq_prune(shared_queue, best_score);
// 			REPORT_ERROR(MPI_Isend(&recv_bound, 1, bound_comm_t, getRightRank(my_rank, torus),
// 									BOUND_TAG, torus, &request),
// 						"MPI_Isend : ERROR\n", );
// 		}
// 	}
//
// 	recv_flag = 0;
// 	REPORT_ERROR(MPI_Iprobe(getRightRank(my_rank, torus), BOUND_TAG, torus, &recv_flag,
// 							&status), "MPI_Iprobe : ERROR\n", );
// 	if (recv_flag) {
// 		REPORT_ERROR(MPI_Recv(&recv_bound, 1, bound_comm_t, getRightRank(my_rank, torus),
// 							  BOUND_TAG, torus, MPI_STATUS_IGNORE),
// 			   "MPI_Recv : ERROR\n", );
//
// 		flag = 0;
// 		omp_set_lock(&best_solution_lock);
// 		if (recv_bound->new_bound < best_score) {
// 			best_score = recv_bound->new_bound;
// 			best_solution = NULL;
// 			flag = 1;
// 		}
// 		omp_unset_lock(&best_solution_lock);
//
// 		if (flag == 1) {
// 			pq_prune(shared_queue, best_score);
// 			REPORT_ERROR(MPI_Isend(&recv_bound, 1, bound_comm_t, getLeftRank(my_rank, torus),
// 								   BOUND_TAG, torus, &request),
// 				"MPI_Isend : ERROR\n", );
// 		}
// 	}
//
// 	recv_flag = 0;
// 	REPORT_ERROR(MPI_Iprobe(getDownRank(my_rank, torus), BOUND_TAG, torus, &recv_flag,
// 							&status), "MPI_Iprobe : ERROR\n", );
// 	if (recv_flag) {
// 		REPORT_ERROR(MPI_Recv(&recv_bound, 1, bound_comm_t, getDownRank(my_rank, torus),
// 							  BOUND_TAG, torus, MPI_STATUS_IGNORE),
// 			   "MPI_Recv : ERROR\n", );
//
// 		flag = 0;
// 		omp_set_lock(&best_solution_lock);
// 		if (recv_bound->new_bound < best_score) {
// 			best_score = recv_bound->new_bound;
// 			best_solution = NULL;
// 			flag = 1;
// 		}
// 		omp_unset_lock(&best_solution_lock);
//
// 		if (flag == 1) {
// 			pq_prune(shared_queue, best_score);
// 			REPORT_ERROR(MPI_Isend(&recv_bound, 1, bound_comm_t, getUpRank(my_rank, torus),
// 								   BOUND_TAG, torus, &request2),
// 				"MPI_Isend : ERROR\n", );
// 			REPORT_ERROR(MPI_Isend(&recv_bound, 1, bound_comm_t, getLeftRank(my_rank, torus),
// 								   BOUND_TAG, torus, &request),
// 				"MPI_Isend : ERROR\n", );
// 			REPORT_ERROR(MPI_Isend(&recv_bound, 1, bound_comm_t, getRightRank(my_rank, torus),
// 								   BOUND_TAG, torus, &request),
// 				"MPI_Isend : ERROR\n", );
// 		}
// 	}
//
// 	recv_flag = 0;
// 	REPORT_ERROR(MPI_Iprobe(getUpRank(my_rank, torus), BOUND_TAG, torus, &recv_flag,
// 							&status), "MPI_Iprobe : ERROR\n", );
// 	if (recv_flag) {
// 		REPORT_ERROR(MPI_Recv(&recv_bound, 1, bound_comm_t, getUpRank(my_rank, torus),
// 							  BOUND_TAG, torus, MPI_STATUS_IGNORE),
// 			   "MPI_Recv : ERROR\n", );
//
// 		flag = 0;
// 		omp_set_lock(&best_solution_lock);
// 		if (recv_bound->new_bound < best_score) {
// 			best_score = recv_bound->new_bound;
// 			best_solution = NULL;
// 			flag = 1;
// 		}
// 		omp_unset_lock(&best_solution_lock);
//
// 		if (flag == 1) {
// 			pq_prune(shared_queue, best_score);
// 			REPORT_ERROR(MPI_Isend(&recv_bound, 1, bound_comm_t, getDownRank(my_rank, torus),
// 								   BOUND_TAG, torus, &request2),
// 				"MPI_Isend : ERROR\n", );
// 			REPORT_ERROR(MPI_Isend(&recv_bound, 1, bound_comm_t, getLeftRank(my_rank, torus),
// 								   BOUND_TAG, torus, &request),
// 				"MPI_Isend : ERROR\n", );
// 			REPORT_ERROR(MPI_Isend(&recv_bound, 1, bound_comm_t, getRightRank(my_rank, torus),
// 								   BOUND_TAG, torus, &request),
// 				"MPI_Isend : ERROR\n", );
// 		}
// 	}
// }

// not threadsafe. use with caution.
void pack_array(queue_head* qh, void* outbuff, int buff_size, int* pos, void* problem_data) {
	if (qh != NULL) {
		pack_solution(outbuff, buff_size, pos, qh->partial_solution, qh->priority, torus,
					  problem_data);
		pack_array(qh->left_subtree, outbuff, buff_size, pos, problem_data);
		pack_array(qh->right_subtree, outbuff, buff_size, pos, problem_data);
	}
}

// not threadsafe. use with caution.
void unpack_array(queue* q, void* outbuff, int buff_size, int len, int* pos,
				  void* problem_data) {
	float score;
	solution_vector psol;

	for (int i = 0; i < len; i++) {
		psol = unpack_solution(outbuff, buff_size, torus, pos, &score, problem_data);
		pq_insert_nc(q, score, psol);
	}

}

int initializeLoad(void* domain_data) {
	if (my_rank == 0) {
		pq_insert_nc(shared_queue, 0.0,
					 get_root_partial_solution(domain_data));

		#pragma omp parallel num_threads(NUM_THREADS)
		{
			queue* priv_queue = create_queue();
			while (pq_length(shared_queue) < X_LIM*Y_LIM) {
				expand_partial_solution(priv_queue, domain_data);

				#pragma omp barrier
				if (pq_length(shared_queue) == 0) {
					end_program = 1;
					break;
				}
			}
		}

		if (end_program) {
			for (int i = 1; i < world_size; i++) {
				MPI_Send(&end_program, 1, MPI_INT, i, INIT_LOAD_TAG, torus);
			}
			// cleanup and shutdown
			return 1;
		} else {
			queue_head* qh;
			int pos;
			void* outbuff = malloc(solution_vector_size);

			for (int i = 1; i < world_size; i++) {
				qh = pq_extract(shared_queue, 1);
				pack_array(qh, outbuff, solution_vector_size, &pos, domain_data);

				MPI_Send(&end_program, 1, MPI_INT, i, INIT_LOAD_TAG, torus);
				MPI_Send(outbuff, solution_vector_size, MPI_PACKED, i, INIT_LOAD_TAG, torus);
			}

			free(outbuff);

			return 0;
		}
	} else {
		int len;

		MPI_Recv(&end_program, 1, MPI_INT, 0, INIT_LOAD_TAG, torus, MPI_STATUS_IGNORE);

		if (end_program) {
			return 1;
		} else {
			float score;
			int pos = 0;
			void* outbuff = malloc(solution_vector_size);
			MPI_Recv(outbuff, solution_vector_size, MPI_PACKED, 0, INIT_LOAD_TAG, torus,
					 MPI_STATUS_IGNORE);
			unpack_array(shared_queue, outbuff, solution_vector_size, 1, &pos,
						 domain_data);

			free(outbuff);

			return 0;
		}
	}
}
