#include <stdlib.h>
#include <stdio.h>

#include "list.h"
#include "comforts.h"
#include "generics.h"

list* create_list() {
	NEW(list, l);
	l->len = 0;
	l->head = NULL;

	return l;
}

void list_insert(list* l, void* item) {
	NEW(struct list_node, le);
	le->val = item;
	le->next = l->head;
	l->head = le;
	l->len++;
}

void list_clear(list* l) {
	l->len = 0;

	struct list_node* le = l->head;
	struct list_node* lf;

	while (le != 0) {
		lf = le->next;
		free(le->val);
		free(le);
		le = lf;
	}
	l->head = NULL;
}

void list_print(list* l, float score) {
	struct list_node* le = l->head;
	struct list_node* lf;

	while (le != 0) {
		print_solution(le->val, score);
		printf("\n");
		lf = le->next;
		free(le->val);
		free(le);
		le = lf;
	}
}
