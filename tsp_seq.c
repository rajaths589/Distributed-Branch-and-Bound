#include "graph.h"
#include "bitvector.h"
#include "priority_queue.h"

#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <float.h>
#include <stddef.h>

#define NUM_THREADS 4

int max_length;
int solution_vector_size = 0;

int working[NUM_THREADS];
int can_work[NUM_THREADS];
omp_lock_t work_lock_turnstile;
queue* shared_queue;

omp_lock_t best_solution_lock;
int best_score_changed = 0;
solution_vector best_solution;
float best_score = FLT_MAX;

/* from generics.h */

typedef void* solution_vector;
extern int solution_vector_size;

void* populate_domain_data(int argc, char** argv);

solution_vector get_root_partial_solution(void* domain_specific_data);

// returns if the partial_solution was a full solution
// i.e., if 1 is returned, partial_solution is a valid solution_vector
// else extension is possible.
int construct_candidates(solution_vector partial_solution, float partial_soln_score,
                                                                queue* private_queue, void* domain_specific_data);

/* from comforts.h */
void print_solution(solution_vector solution, float score);

#define NEW(TYPE, X) TYPE* X = (TYPE*) malloc(sizeof(TYPE))
#define MIN(X, Y) X<Y ? X:Y
#define MAX(X, Y) X>Y ? X:Y

// assumes that temp is declared
#define SWAP(X, Y) temp = X; X = Y; Y = temp

extern int err;

typedef struct tsp_data {
        graph* g;
        int start_point;
        int max_len;	// equal to num vertices
} tsp_data;

struct bound_array {
        int* path;
        int curr_length;
        int max_length;
        bitvector* used_vertices;
};

typedef struct bound_array tsp_path;

tsp_path* create_soln_copy(tsp_path* p) {
        NEW(tsp_path, p1);
        p1->curr_length = p->curr_length;
        p1->max_length = p->max_length;
        p1->path = calloc(p1->max_length, sizeof(int));
        memcpy(p1->path, p->path, p1->max_length*sizeof(int));
        p1->used_vertices = create_copy(p->used_vertices);

        return p1;
}

void* populate_domain_data(int argc, char** argv) {
        assert(argc == 2);
        FILE* fp = fopen(argv[1], "r");

        assert(fp != NULL);

        int num_vertices, num_edges;
        NEW(tsp_data, data);

        fscanf(fp, "%d %d", &num_vertices, &num_edges);
        graph* g = create_graph(num_vertices);

        int from, to;
        float weight;

        for (int i = 0; i < num_edges; i++) {
                fscanf(fp, "%d %d %f", &from, &to, &weight);
                add_edge(g, from, to, weight);
        }

        data->g = g;
        data->max_len = g->num_vertices;

        fscanf(fp, "%d", &data->start_point);

        assert(data->start_point >= 0);
        assert(data->start_point < num_vertices);

        solution_vector_size = sizeof(float) + sizeof(int) + (data->max_len+1)*sizeof(int);

        return (void*) data;
}

solution_vector get_root_partial_solution(void* domain_specific_data) {
        NEW(tsp_path, empty_solution);

        tsp_data* data = (tsp_data*) domain_specific_data;

        empty_solution->max_length = data->max_len + 1;
        empty_solution->curr_length = 1;
        empty_solution->path = (int*) calloc(data->max_len + 1, sizeof(int));
        empty_solution->path[0] = data->start_point;
        empty_solution->used_vertices = create_bitvector(data->max_len);
        setIndex(empty_solution->used_vertices, data->start_point);

        return (void*) empty_solution;
}

int construct_candidates(solution_vector partial_solution, float partial_soln_score,
                                                 queue* private_queue, void* domain_specific_data) {

        tsp_path* partial = (tsp_path*) partial_solution;
        tsp_data* data = (tsp_data*) domain_specific_data;

        if (partial->curr_length == data->max_len+1)
                return 1;

        list_node* l;
        tsp_path* path1;
        l = data->g->adjacency_list[partial->path[partial->curr_length-1]];

        if (partial->curr_length == data->max_len) {
                while (l != NULL) {
                        if (l->to == data->start_point) {
                                partial->path[partial->curr_length] = l->to;
                                partial->curr_length++;
                                pq_insert_nc(private_queue, partial_soln_score + l->weight, partial);
                                break;
                        }
                        l = l->next;
                }
        } else {
                while (l != NULL) {
                        if (!getIndex(partial->used_vertices, l->to)) {
                                path1 = create_soln_copy(partial);
                                path1->path[path1->curr_length] = l->to;
                                path1->curr_length++;
                                setIndex(path1->used_vertices, l->to);
                                pq_insert_nc(private_queue, partial_soln_score + l->weight, path1);
                        }
                        l = l->next;
                }

                free(partial);
        }

        return 0;
}

void print_solution(solution_vector sol, float score) {
        tsp_path* solution = (tsp_path*) sol;

        printf("Curr length : %d\n", solution->curr_length);
        for (int i = 0; i < solution->curr_length; i++) {
                printf("%d\t", solution->path[i]);
        }
        printf("\n");
        printf("Cost : %f\n", score);
}

int safe_solution_buffer_size(tsp_data* data) {
        return (data->max_len+1)*sizeof(int) + sizeof(float) + sizeof(int);
}

void expand_partial_solution(queue* private_queue, void* domain_data) {
        float score;

        solution_vector partial_solution = pq_min_extract(shared_queue, &score);

        if (partial_solution != NULL) {
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
        }
}

int main(int argc, char** argv) {
        omp_init_lock(&best_solution_lock);
        omp_init_lock(&work_lock_turnstile);

        shared_queue = create_queue();
        for (int i = 0; i < NUM_THREADS; i++) {
                working[i] = 0;
                can_work[i] = 0;
        }

        void* domain_data = populate_domain_data(argc, argv);

        pq_insert_nc(shared_queue, 0.0,
                get_root_partial_solution(domain_data));

        #pragma omp parallel num_threads(NUM_THREADS)
        {
                queue* priv_queue = create_queue();
                int l;
                while (1) {
                        l = pq_length(shared_queue);
                        #pragma omp barrier
                        if (l == 0) {
                                break;
                        }

                        expand_partial_solution(priv_queue, domain_data);
                        #pragma omp barrier
                }

                destroy_queue(priv_queue);
        }

        #pragma omp parallel num_threads(NUM_THREADS)
        {
                int thread_rank = omp_get_thread_num();

                if (thread_rank != NUM_THREADS-1) {
                        queue* private_queue = create_queue();
                        int l;
                        do {
                                omp_set_lock(&work_lock_turnstile);
                                omp_unset_lock(&work_lock_turnstile);

                                working[thread_rank] = 1;

                                expand_partial_solution(private_queue, domain_data);

                                working[thread_rank] = 0;
                                l = pq_length(private_queue);
                        } while (l);

                }
        }
}
