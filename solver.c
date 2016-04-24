#include <omp.h>
#include <mpi.h>
#include <float.h>

#include "priority_queue.h"
#include "generics.h"
#include "comforts.h"

typedef struct bound_comm {
	float new_bound;
	int start_rank;
} bound_comm;

MPI_Datatype bound_comm_t;

omp_lock_t best_solution_lock;
int best_score_changed = 0;
float best_score = FLT_MAX;
solution_vector best_solution;

void expand_partial_solution(queue* shared_queue, queue* private_queue, void* domain_data) {
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

void send_bound() {
	MPI_Request request;
	bound_comm recv_bound;
	recv_bound->start_rank = my_rank;

	omp_set_lock(&best_solution_lock);
	recv_bound->new_bound = best_score;
	best_score_changed = 0;
	omp_unset_lock(&best_solution_lock);

	REPORT_ERROR(MPI_Isend(&recv_bound, 1, bound_comm_t,  getUpRank(my_rank, torus),
				BOUND_TAG, torus, &request), "MPI_Isend : ERROR\n", );
	REPORT_ERROR(MPI_Isend(&recv_bound, 1, bound_comm_t,  getDownRank(my_rank, torus),
						   BOUND_TAG, torus, &request), "MPI_Isend : ERROR\n", );
	REPORT_ERROR(MPI_Isend(&recv_bound, 1, bound_comm_t,  getLeftRank(my_rank, torus),
						   BOUND_TAG, torus, &request), "MPI_Isend : ERROR\n", );
	REPORT_ERROR(MPI_Isend(&recv_bound, 1, bound_comm_t,  getRightRank(my_rank, torus),
						   BOUND_TAG, torus, &request), "MPI_Isend : ERROR\n", );
}

void receive_and_forward_bound() {
	int recv_flag = 0;
	int flag = 0;
	MPI_Status status;
	MPI_Request request;
	bound_comm recv_bound;

	REPORT_ERROR(MPI_Iprobe(getLeftRank(my_rank, torus), BOUND_TAG, torus, &recv_flag,
							&status), "MPI_Iprobe : ERROR\n", );
	if (recv_flag) {
		REPORT_ERROR(MPI_Recv(&recv_bound, 1, bound_comm_t, getLeftRank(my_rank, torus),
							BOUND_TAG, torus, MPI_STATUS_IGNORE),
					"MPI_Recv : ERROR\n", );

		flag = 0;
		omp_set_lock(&best_solution_lock);
		if (recv_bound->new_bound < best_score) {
			best_score = recv_bound->new_bound;
			best_solution = NULL;
			flag = 1;
		}
		omp_unset_lock(&best_solution_lock);

		if (flag == 1) {
			pq_prune(shared_queue, best_score);
			REPORT_ERROR(MPI_Isend(&recv_bound, 1, bound_comm_t, getRightRank(my_rank, torus),
									BOUND_TAG, torus, &request),
						"MPI_Isend : ERROR\n", );
		}
	}

	recv_flag = 0;
	REPORT_ERROR(MPI_Iprobe(getRightRank(my_rank, torus), BOUND_TAG, torus, &recv_flag,
							&status), "MPI_Iprobe : ERROR\n", );
	if (recv_flag) {
		REPORT_ERROR(MPI_Recv(&recv_bound, 1, bound_comm_t, getRightRank(my_rank, torus),
							  BOUND_TAG, torus, MPI_STATUS_IGNORE),
			   "MPI_Recv : ERROR\n", );

		flag = 0;
		omp_set_lock(&best_solution_lock);
		if (recv_bound->new_bound < best_score) {
			best_score = recv_bound->new_bound;
			best_solution = NULL;
			flag = 1;
		}
		omp_unset_lock(&best_solution_lock);

		if (flag == 1) {
			pq_prune(shared_queue, best_score);
			REPORT_ERROR(MPI_Isend(&recv_bound, 1, bound_comm_t, getLeftRank(my_rank, torus),
								   BOUND_TAG, torus, &request),
				"MPI_Isend : ERROR\n", );
		}
	}

	recv_flag = 0;
	REPORT_ERROR(MPI_Iprobe(getDownRank(my_rank, torus), BOUND_TAG, torus, &recv_flag,
							&status), "MPI_Iprobe : ERROR\n", );
	if (recv_flag) {
		REPORT_ERROR(MPI_Recv(&recv_bound, 1, bound_comm_t, getDownRank(my_rank, torus),
							  BOUND_TAG, torus, MPI_STATUS_IGNORE),
			   "MPI_Recv : ERROR\n", );

		flag = 0;
		omp_set_lock(&best_solution_lock);
		if (recv_bound->new_bound < best_score) {
			best_score = recv_bound->new_bound;
			best_solution = NULL;
			flag = 1;
		}
		omp_unset_lock(&best_solution_lock);

		if (flag == 1) {
			pq_prune(shared_queue, best_score);
			REPORT_ERROR(MPI_Isend(&recv_bound, 1, bound_comm_t, getUpRank(my_rank, torus),
								   BOUND_TAG, torus, &request2),
				"MPI_Isend : ERROR\n", );
			REPORT_ERROR(MPI_Isend(&recv_bound, 1, bound_comm_t, getLeftRank(my_rank, torus),
								   BOUND_TAG, torus, &request),
				"MPI_Isend : ERROR\n", );
			REPORT_ERROR(MPI_Isend(&recv_bound, 1, bound_comm_t, getRightRank(my_rank, torus),
								   BOUND_TAG, torus, &request),
				"MPI_Isend : ERROR\n", );
		}
	}

	recv_flag = 0;
	REPORT_ERROR(MPI_Iprobe(getUpRank(my_rank, torus), BOUND_TAG, torus, &recv_flag,
							&status), "MPI_Iprobe : ERROR\n", );
	if (recv_flag) {
		REPORT_ERROR(MPI_Recv(&recv_bound, 1, bound_comm_t, getUpRank(my_rank, torus),
							  BOUND_TAG, torus, MPI_STATUS_IGNORE),
			   "MPI_Recv : ERROR\n", );

		flag = 0;
		omp_set_lock(&best_solution_lock);
		if (recv_bound->new_bound < best_score) {
			best_score = recv_bound->new_bound;
			best_solution = NULL;
			flag = 1;
		}
		omp_unset_lock(&best_solution_lock);

		if (flag == 1) {
			pq_prune(shared_queue, best_score);
			REPORT_ERROR(MPI_Isend(&recv_bound, 1, bound_comm_t, getDownRank(my_rank, torus),
								   BOUND_TAG, torus, &request2),
				"MPI_Isend : ERROR\n", );
			REPORT_ERROR(MPI_Isend(&recv_bound, 1, bound_comm_t, getLeftRank(my_rank, torus),
								   BOUND_TAG, torus, &request),
				"MPI_Isend : ERROR\n", );
			REPORT_ERROR(MPI_Isend(&recv_bound, 1, bound_comm_t, getRightRank(my_rank, torus),
								   BOUND_TAG, torus, &request),
				"MPI_Isend : ERROR\n", );
		}
	}
}


void queue_balancing(int *loads_arr, MPI_Comm current_comm );

void load_balancing(MPI_Comm *load_comm){
	int stop_flag, share, i, queue_len, flag[5], loads_arr[5];
	share = 1;

	MPI_Request req_comm[5];
	MPI_Status status_comm[5];

	for(i=0; i<5; i++){
		//i=0: I am the Gatherer
		REPORT_ERROR(MPI_Ibcast(&share, 1, MPI_INT, CENTRE, load_comm[i], &req_comm[i]),
					"MPI_Ibcast : ERROR\n", );
	}

	int order[4] = {1,2,3,4};

	while(stop_flag){

		i = 0;
		REPORT_ERROR(MPI_Test(&req_comm[CENTRE], &flag[CENTRE], &status_comm[CENTRE]),
					"MPI_Test[CENTRE] : ERROR\n", );

		if(flag[CENTRE]){

			queue_len = getQueueLen();

			REPORT_ERROR(MPI_Gather(&queue_len, 1, MPI_INT, loads_arr, 1, MPI_INT, CENTRE, load_comm[order[i]]),
				"MPI_Gather : ERROR\n", );

			queue_balancing(loads_arr, load_comm[CENTRE]);

			//EVENTUALLY:
			REPORT_ERROR(MPI_Ibcast(&share, 1, MPI_INT, CENTRE, load_comm[CENTRE], &req_comm[CENTRE]),
				"MPI_Ibcast : ERROR\n", );

		}

		random_perm(order, 4);

		for(i=0; i<4; i++){

			REPORT_ERROR(MPI_Test(&req_comm[order[i]], &flag[order[i]], &status_comm[order[i]]),
					"MPI_Test : ERROR\n", );

			if(flag[order[i]]){

				//DO SOMETHING:
				queue_len = getQueueLen();
				REPORT_ERROR(MPI_Gather(&queue_len, 1, MPI_INT, NULL, 0, MPI_INT, CENTRE, load_comm[order[i]]),
					"MPI_Gather : ERROR\n", );

				queue_balancing(loads_arr, load_comm[order[i]]);

				//EVENTUALLY:
				REPORT_ERROR(MPI_Ibcast(&share, 1, MPI_INT, CENTRE, load_comm[i], &req_comm[i]),
					"MPI_Ibcast : ERROR\n", );

			}
		}
	}
}


void queue_balancing(int *loads_arr, MPI_Comm current_comm ){

	int my_rank;
	MPI_Comm_rank(current_comm, &my_rank);

	if(my_rank == CENTRE){

	}

	else {

	}
}