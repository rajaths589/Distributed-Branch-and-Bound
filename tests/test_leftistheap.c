#include <omp.h>
#include "priority_queue.h"
#include <stdlib.h>
#include <time.h>
#include <stdio.h>

int main(int argc, char* argv[]) {
	int nthreads = 4;
	queue* q = create_queue();

	srand(time(NULL));
	int max_numbers = 1000;
	int num_inserts = 100;

//	unsigned char* is_used = (unsigned char*) calloc(max_numbers, sizeof(unsigned char));
	int* numbers = (int*) malloc(num_inserts*sizeof(int));

	for (int i = 0; i < num_inserts; ) {
		int r = rand()%max_numbers;

//		if (is_used[r])
//			continue;

//		is_used[r] = 1;
		numbers[i] = r;
		i++;
	}

	int p = num_inserts/nthreads;

	#pragma omp parallel num_threads(nthreads)
	{
		int my_rank = omp_get_thread_num();

		queue* private_q = create_queue();

		for (int j = my_rank*p; j < (my_rank+1)*p; j++) {
			pq_insert_nc(private_q, numbers[j], NULL);
		}

		pq_merge(q, private_q);
	}

	pq_prune(q, 200);

	float f;
	while (pq_length(q)) {
		pq_min_extract(q, &f);

		printf("%f\n", f);
	}

	return 0;
}
