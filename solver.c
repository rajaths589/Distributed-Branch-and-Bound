#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
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
#define BOUND_TAG 55

#define NUM_THREADS 3

#define MIN_NULL_FLAG_LOAD_BAL 5
#define MAX_QLEN_LOAD_BAL 80

int working[NUM_THREADS];
int null_flag;
omp_lock_t work_lock_turnstile;
queue* shared_queue;

int err = 0;
int my_rank;
int world_size;
int X_LIM;
int Y_LIM;
int end_program = 0;
int in_loadbal = 0;

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

void send_bound();
void receive_and_forward_bound();
void queue_balancing(int neighbor_rank, int my_queue_len, int neighbour_queue_len);
void expand_partial_solution(queue* private_queue, void* domain_data);
MPI_Comm* create_topology();
void setupCommunicators();
void load_balancing();
void loadbal_recipient(enum CommNum direc);
void loadbal_initiator(enum CommNum direc);
int initializeLoad(void* domain_data);

MPI_Datatype bound_comm_t;

omp_lock_t best_solution_lock;
int best_score_changed = 0;
float best_score = FLT_MAX;
solution_vector best_solution;
MPI_Comm torus;
int neighbors_rank[5];
int torus_neighbors_rank[5];
int stop_flag = 0;


int main(int argc, char** argv) {
	MPI_Init(&argc, &argv);

	setupCommunicators();
	MPI_Comm_rank(torus, &my_rank);
	srand(time(NULL));

	const int dt_num = 2;
	MPI_Datatype dt_type[2] = {MPI_FLOAT, MPI_INT};
	int dt_blocklen[2] = {1, 1};
	MPI_Aint offset[2];
	offset[0] = offsetof(bound_comm, new_bound);
	offset[1] = offsetof(bound_comm, start_rank);

	MPI_Type_create_struct(dt_num, dt_blocklen, offset, dt_type, &bound_comm_t);
	MPI_Type_commit(&bound_comm_t);

	omp_init_lock(&best_solution_lock);
	omp_init_lock(&work_lock_turnstile);

	void* domain_data = populate_domain_data(argc, argv);

	shared_queue = create_queue();

	initializeLoad(domain_data);
	null_flag = 0;

	#pragma omp parallel num_threads(NUM_THREADS)
	{
		int thread_rank = omp_get_thread_num();

		if (thread_rank != NUM_THREADS-1) {
			queue* private_queue = create_queue();
			working[thread_rank] = 0;

			while (!end_program) {
				omp_set_lock(&work_lock_turnstile);
				omp_unset_lock(&work_lock_turnstile);

				#pragma omp atomic
				working[thread_rank] |= 1;

//				printf("RANK:%d THREAD:%d\n", my_rank, thread_rank);
				expand_partial_solution(private_queue, domain_data);

				#pragma omp atomic
				working[thread_rank] &= 0;
			}

		} else {
			int bound_send_flag, flag;
			int syn_receive, syn_ack_receive, syn_share=1;
			MPI_Request syn_request_down, syn_request_right;
			MPI_Request syn_receive_left, syn_receive_up, syn_ack_right, syn_ack_down;

			if (neighbors_rank[LEFT] != -1) {
				MPI_Irecv(&syn_receive, 1, MPI_INT, neighbors_rank[LEFT], SYN_TAG, torus,
							&syn_receive_left);
			}

			if (neighbors_rank[UP] != -1) {
				MPI_Irecv(&syn_receive, 1, MPI_INT, neighbors_rank[UP], SYN_TAG, torus,
							&syn_receive_up);
			}

			while (!end_program) {
				bound_send_flag = 0;
				omp_set_lock(&best_solution_lock);
				if (best_score_changed > 0) {
					bound_send_flag = 1;
				}
				omp_unset_lock(&best_solution_lock);

				if (bound_send_flag)
					send_bound();

				receive_and_forward_bound();

				if (neighbors_rank[LEFT] != -1) {
					MPI_Test(&syn_receive_left, &flag, MPI_STATUS_IGNORE);

					if (flag) {
						loadbal_recipient(LEFT);
						MPI_Irecv(&syn_receive, 1, MPI_INT, neighbors_rank[LEFT], SYN_TAG,
								  torus, &syn_receive_left);
					}
				}

				if (neighbors_rank[UP] != -1) {
					MPI_Test(&syn_receive_up, &flag, MPI_STATUS_IGNORE);

					if (flag) {
						loadbal_recipient(UP);
						MPI_Irecv(&syn_receive, 1, MPI_INT, neighbors_rank[UP], SYN_TAG,
								torus, &syn_receive_up);
					}
				}

				if (in_loadbal) {
					if (neighbors_rank[RIGHT] != -1) {
						MPI_Test(&syn_ack_right, &flag, MPI_STATUS_IGNORE);

						if (flag) {
							loadbal_initiator(RIGHT);
							in_loadbal--;
						}
					}

					if (neighbors_rank[DOWN] != -1) {
						MPI_Test(&syn_ack_down, &flag, MPI_STATUS_IGNORE);

						if (flag) {
							loadbal_initiator(DOWN);
							in_loadbal--;
						}
					}
				}

				if (!in_loadbal && null_flag >= MIN_NULL_FLAG_LOAD_BAL) {
					//do load_balancing

					if (neighbors_rank[RIGHT] != -1) {
						in_loadbal++;
						MPI_Isend(&syn_share, 1, MPI_INT, neighbors_rank[RIGHT], SYN_TAG,
								  torus, &syn_request_right);
						MPI_Irecv(&syn_ack_receive, 1, MPI_INT,
								  neighbors_rank[RIGHT],SYN_ACK_TAG,torus, &syn_ack_right);
					}
					if (neighbors_rank[DOWN] != -1) {
						in_loadbal++;
						MPI_Isend(&syn_share, 1, MPI_INT, neighbors_rank[DOWN], SYN_TAG,
								  torus, &syn_request_down);
						MPI_Irecv(&syn_ack_receive, 1, MPI_INT,
								  neighbors_rank[DOWN],SYN_ACK_TAG,torus, &syn_ack_down);
					}

					#pragma omp atomic
					null_flag ^= null_flag;
				} else if (!in_loadbal && pq_length(shared_queue) > MAX_QLEN_LOAD_BAL) {
					//do load_balancing

					if (neighbors_rank[RIGHT] != -1) {
						in_loadbal++;
						MPI_Isend(&syn_share, 1, MPI_INT, neighbors_rank[RIGHT], SYN_TAG,
								  torus, &syn_request_right);
						MPI_Irecv(&syn_ack_receive, 1, MPI_INT,
									neighbors_rank[RIGHT],SYN_ACK_TAG,torus, &syn_ack_right);
					}
					if (neighbors_rank[DOWN] != -1) {
						in_loadbal++;
						MPI_Isend(&syn_share, 1, MPI_INT, neighbors_rank[DOWN], SYN_TAG,
								  torus, &syn_request_down);
						MPI_Irecv(&syn_ack_receive, 1, MPI_INT,
								  neighbors_rank[DOWN],SYN_ACK_TAG,torus, &syn_ack_down);
					}

					#pragma omp atomic
					null_flag ^= null_flag;
				}
			}
		}
	}

	MPI_Finalize();

	return 0;
}

void loadbal_recipient(enum CommNum direc) {
	int syn_ack_send = 1;
	int send_queue_len, recv_queue_len;

	MPI_Send(&syn_ack_send, 1, MPI_INT, neighbors_rank[direc],
			 SYN_ACK_TAG, torus);
	MPI_Recv(&recv_queue_len, 1, MPI_INT, neighbors_rank[direc], DATA_TAG,
			 torus, MPI_STATUS_IGNORE);

	omp_set_lock(&work_lock_turnstile);
	int all_done = 0;
	while (!all_done) {
		for (int i = 0; i < NUM_THREADS-1; i++) {
			#pragma omp atomic
			all_done |= working[i];

			if (all_done)
				break;
		}

		if (all_done)
			all_done = 0;
		else
			break;
	}
	send_queue_len = pq_length(shared_queue);
	omp_unset_lock(&work_lock_turnstile);

	MPI_Send(&send_queue_len, 1, MPI_INT, neighbors_rank[direc],DATA_TAG,
			 torus);
	queue_balancing(neighbors_rank[direc], send_queue_len, recv_queue_len);
}

void loadbal_initiator(enum CommNum direc) {
	int send_queue_len, recv_queue_len;

	omp_set_lock(&work_lock_turnstile);
	int all_done = 0;
	while (!all_done) {
		for (int i = 0; i < NUM_THREADS-1; i++) {
			#pragma omp atomic
			all_done |= working[i];

			if (all_done)
				break;
		}

		if (all_done)
			all_done = 0;
		else
			break;
	}
	send_queue_len = pq_length(shared_queue);
	omp_unset_lock(&work_lock_turnstile);

	MPI_Send(&send_queue_len, 1, MPI_INT, neighbors_rank[direc], DATA_TAG, torus);
	MPI_Recv(&recv_queue_len, 1, MPI_INT, neighbors_rank[direc], DATA_TAG, torus,
			MPI_STATUS_IGNORE);

	queue_balancing(neighbors_rank[direc], send_queue_len, recv_queue_len);
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
	int coords[2];

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

		MPI_Cart_rank(torus, coords_ptr[i], &torus_neighbors_rank[i]);
	}
}

MPI_Comm* create_topology() {
	MPI_Comm *mesh = (MPI_Comm *) malloc(sizeof(MPI_Comm));
	int ndims = 2;
	int dim_sz[2];
	int periodic[2];

	dim_sz[0] = X_LIM;
	dim_sz[1] = Y_LIM;

	periodic[0] = 1;
	periodic[1] = 1;

	int reorder = 1;

	REPORT_ERROR(MPI_Cart_create(MPI_COMM_WORLD, ndims, dim_sz, periodic, reorder, mesh),
				 "MPI_Cart_create: ERROR\n",
			  NULL);

	return mesh;
}

void expand_partial_solution(queue* private_queue, void* domain_data) {
	float score;

	solution_vector partial_solution = pq_min_extract(shared_queue, &score);

	if (partial_solution != NULL) {
// 		printf("RANK: %d\n", my_rank);
// 		print_solution(partial_solution, score);
// 		printf("\n");

		if (construct_candidates(partial_solution, score, private_queue, domain_data)) {
			int flag = 0;

			omp_set_lock(&best_solution_lock);
			if (score < best_score) {
				best_score = score;
				best_solution = partial_solution;
				printf("BEST: \n");
				print_solution(partial_solution, best_score);
				printf("\n");
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
	} else {
		#pragma omp atomic
		null_flag += 1;
	}
}

// int getQueueLen() {
// 	return rand()%100;
// }

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
void queue_balancing(int n_rank, int my_queue_len, int neighbour_queue_len) {
	printf("RANK:%d NEIGHBOUR_RANK:%d MY_QUEUE:%d NEIGHBOUR_QUEUE:%d \n", my_rank,
		   n_rank, my_queue_len, neighbour_queue_len);

}

void send_bound() {
	MPI_Request request;
	bound_comm recv_bound;
	recv_bound.start_rank = my_rank;

	omp_set_lock(&best_solution_lock);
	recv_bound.new_bound = best_score;
	best_score_changed = 0;
	omp_unset_lock(&best_solution_lock);

	REPORT_ERROR(MPI_Isend(&recv_bound, 1, bound_comm_t,  torus_neighbors_rank[UP],
				BOUND_TAG, torus, &request), "MPI_Isend : ERROR\n", );
	REPORT_ERROR(MPI_Isend(&recv_bound, 1, bound_comm_t,  torus_neighbors_rank[DOWN],
						   BOUND_TAG, torus, &request), "MPI_Isend : ERROR\n", );
	REPORT_ERROR(MPI_Isend(&recv_bound, 1, bound_comm_t,  torus_neighbors_rank[LEFT],
						   BOUND_TAG, torus, &request), "MPI_Isend : ERROR\n", );
	REPORT_ERROR(MPI_Isend(&recv_bound, 1, bound_comm_t,  torus_neighbors_rank[RIGHT],
						   BOUND_TAG, torus, &request), "MPI_Isend : ERROR\n", );
}

void receive_and_forward_bound() {
	int recv_flag = 0;
	int flag = 0;
	MPI_Status status;
	MPI_Request request, request2;
	bound_comm recv_bound;

	REPORT_ERROR(MPI_Iprobe(torus_neighbors_rank[LEFT], BOUND_TAG, torus, &recv_flag,
							&status), "MPI_Iprobe : ERROR\n", );
	if (recv_flag) {
		REPORT_ERROR(MPI_Recv(&recv_bound, 1, bound_comm_t, torus_neighbors_rank[LEFT],
							BOUND_TAG, torus, MPI_STATUS_IGNORE),
					"MPI_Recv : ERROR\n", );
// 		printf("RANK:%d RECVD from LEFT: %f\n", my_rank, recv_bound.new_bound);
		flag = 0;
		omp_set_lock(&best_solution_lock);
		if (recv_bound.new_bound < best_score) {
			best_score = recv_bound.new_bound;
			best_solution = NULL;
			flag = 1;
		}
		omp_unset_lock(&best_solution_lock);

		if (flag == 1) {
			pq_prune(shared_queue, best_score);
			REPORT_ERROR(MPI_Isend(&recv_bound, 1, bound_comm_t, torus_neighbors_rank[RIGHT],
									BOUND_TAG, torus, &request),
						"MPI_Isend : ERROR\n", );
		}
	}

	recv_flag = 0;
	REPORT_ERROR(MPI_Iprobe(torus_neighbors_rank[RIGHT], BOUND_TAG, torus, &recv_flag,
							&status), "MPI_Iprobe : ERROR\n", );
	if (recv_flag) {
		REPORT_ERROR(MPI_Recv(&recv_bound, 1, bound_comm_t, torus_neighbors_rank[RIGHT],
							  BOUND_TAG, torus, MPI_STATUS_IGNORE),
			   "MPI_Recv : ERROR\n", );

// 		printf("RANK:%d RECVD from RIGHT: %f\n", my_rank, recv_bound.new_bound);

		flag = 0;
		omp_set_lock(&best_solution_lock);
		if (recv_bound.new_bound < best_score) {
			best_score = recv_bound.new_bound;
			best_solution = NULL;
			flag = 1;
		}
		omp_unset_lock(&best_solution_lock);

		if (flag == 1) {
			pq_prune(shared_queue, best_score);
			REPORT_ERROR(MPI_Isend(&recv_bound, 1, bound_comm_t, torus_neighbors_rank[LEFT],
								   BOUND_TAG, torus, &request),
				"MPI_Isend : ERROR\n", );
		}
	}

	recv_flag = 0;
	REPORT_ERROR(MPI_Iprobe(torus_neighbors_rank[DOWN], BOUND_TAG, torus, &recv_flag,
							&status), "MPI_Iprobe : ERROR\n", );
	if (recv_flag) {
		REPORT_ERROR(MPI_Recv(&recv_bound, 1, bound_comm_t, torus_neighbors_rank[DOWN],
							  BOUND_TAG, torus, MPI_STATUS_IGNORE),
			   "MPI_Recv : ERROR\n", );

// 		printf("RANK:%d RECVD from DOWN: %f\n", my_rank, recv_bound.new_bound);

		flag = 0;
		omp_set_lock(&best_solution_lock);
		if (recv_bound.new_bound < best_score) {
			best_score = recv_bound.new_bound;
			best_solution = NULL;
			flag = 1;
		}
		omp_unset_lock(&best_solution_lock);

		if (flag == 1) {
			pq_prune(shared_queue, best_score);
			REPORT_ERROR(MPI_Isend(&recv_bound, 1, bound_comm_t, torus_neighbors_rank[UP],
								   BOUND_TAG, torus, &request2),
				"MPI_Isend : ERROR\n", );
			REPORT_ERROR(MPI_Isend(&recv_bound, 1, bound_comm_t, torus_neighbors_rank[LEFT],
								   BOUND_TAG, torus, &request),
				"MPI_Isend : ERROR\n", );
			REPORT_ERROR(MPI_Isend(&recv_bound, 1, bound_comm_t, torus_neighbors_rank[RIGHT],
								   BOUND_TAG, torus, &request),
				"MPI_Isend : ERROR\n", );
		}
	}

	recv_flag = 0;
	REPORT_ERROR(MPI_Iprobe(torus_neighbors_rank[UP], BOUND_TAG, torus, &recv_flag,
							&status), "MPI_Iprobe : ERROR\n", );
	if (recv_flag) {
		REPORT_ERROR(MPI_Recv(&recv_bound, 1, bound_comm_t, torus_neighbors_rank[UP],
							  BOUND_TAG, torus, MPI_STATUS_IGNORE),
			   "MPI_Recv : ERROR\n", );

// 		printf("RANK:%d RECVD from UP: %f\n", my_rank, recv_bound.new_bound);

		flag = 0;
		omp_set_lock(&best_solution_lock);
		if (recv_bound.new_bound < best_score) {
			best_score = recv_bound.new_bound;
			best_solution = NULL;
			flag = 1;
		}
		omp_unset_lock(&best_solution_lock);

		if (flag == 1) {
			pq_prune(shared_queue, best_score);
			REPORT_ERROR(MPI_Isend(&recv_bound, 1, bound_comm_t, torus_neighbors_rank[DOWN],
								   BOUND_TAG, torus, &request2),
				"MPI_Isend : ERROR\n", );
			REPORT_ERROR(MPI_Isend(&recv_bound, 1, bound_comm_t, torus_neighbors_rank[LEFT],
								   BOUND_TAG, torus, &request),
				"MPI_Isend : ERROR\n", );
			REPORT_ERROR(MPI_Isend(&recv_bound, 1, bound_comm_t, torus_neighbors_rank[RIGHT],
								   BOUND_TAG, torus, &request),
				"MPI_Isend : ERROR\n", );
		}
	}
}

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
			int l;
			while (1) {
				l = pq_length(shared_queue);

				#pragma omp barrier

				if (l > X_LIM * Y_LIM) {
					break;
				}

				if (l == 0) {
					end_program = 1;
					break;
				}
				expand_partial_solution(priv_queue, domain_data);

				#pragma omp barrier
			}

			destroy_queue(priv_queue);
		}

		if (end_program) {
			for (int i = 1; i < world_size; i++) {
				MPI_Send(&end_program, 1, MPI_INT, i, INIT_LOAD_TAG, torus);
			}
			// cleanup and shutdown
			return 1;
		} else {
			queue_head* qh;
			int pos = 0;
			void* outbuff = malloc(solution_vector_size);

			for (int i = 1; i < world_size; i++) {
				pos = 0;
				qh = pq_extract(shared_queue, 1);
				pack_array(qh, outbuff, solution_vector_size, &pos, domain_data);

				MPI_Send(&end_program, 1, MPI_INT, i, INIT_LOAD_TAG, torus);
				MPI_Send(outbuff, solution_vector_size, MPI_PACKED, i, INIT_LOAD_TAG, torus);
			}

			free(outbuff);

			return 0;
		}
	} else {
		MPI_Recv(&end_program, 1, MPI_INT, 0, INIT_LOAD_TAG, torus, MPI_STATUS_IGNORE);

		if (end_program) {
			return 1;
		} else {
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
