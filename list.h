/* Solution list. -- not general purpose linked list.
 */

#ifndef LLIST_DEFN
#define LLIST_DEFN

struct list_node {
	struct list_node* next;
	void* val;
};

typedef struct list {
	int len;
	struct list_node* head;
} list;

list* create_list();
void list_insert(list* l, void* item);
void list_clear(list* l);
void list_print(list* l, float score);
#endif
