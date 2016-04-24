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
