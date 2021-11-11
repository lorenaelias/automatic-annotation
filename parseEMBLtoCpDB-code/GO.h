#include <stdio.h>
#define MAX_IDENTIFIER 300

typedef struct GO{
	   char *type;
	   char *value;
       struct GO *next;
} GO_list;

GO_list *head;

void GO_add(char *type, char *value);
void GO_report(FILE *out, char *systematic_id);
int GO_find_id(struct GO *, char**);
