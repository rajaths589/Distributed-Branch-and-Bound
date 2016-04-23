#ifndef GENERIC_TEMPLATE_DEFNS
#define GENERIC_TEMPLATE_DEFNS

typedef void* solution_vector;

void* populate_domain_data(int argc, char** argv);

solution_vector get_root_partial_solution(void* domain_specific_data);

// returns if the partial_solution was a full solution
// i.e., if 1 is returned, partial_solution is a valid solution_vector
// else extension is possible.
int construct_candidates(solution_vector partial_solution, float partial_soln_score,
								queue* private_queue, void* domain_specific_data);


void print_solution(solution_vector solution, float score);

#endif
