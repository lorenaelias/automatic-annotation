#include "GO.h"
#include <malloc.h>
#include <string.h>

void GO_add(char *type, char *value){
    struct GO *item=NULL;

	/* Inserting new item */
	item = (struct GO *)malloc(sizeof(struct GO));
	item->type = (char *)malloc(MAX_IDENTIFIER);
	strcpy(item->type,strdup(type));
	item->value = (char *)malloc(MAX_IDENTIFIER);
	strcpy(item->value,strdup(value));
	if(head!=NULL){
		item->next = head->next;
		head->next = item;
    }else{
		head=item;
		head->next = NULL;
	}
}//function

void GO_report(FILE *out, char *systematic_id){
	struct GO *item=NULL, *clean=NULL;
	item = head;
	char *goid=(char *)malloc(MAX_IDENTIFIER);
	goid[0]='\0';
	
	if (out != NULL && goid != NULL ){
		while (item !=NULL){
			if( !strcmp(item->type,"begin") ){
				GO_find_id(item, &goid);
				item = item->next;
			}
			fprintf(out,"insert into GO (GENE_systematic_id, go_id, go_type, go_value)");
			fprintf(out," values ('%s', '%s', '%s', '%s');\n", systematic_id, goid, item->type, item->value);
			clean = item;
			item = item->next;

			if(clean!=NULL && NULL){
				if(clean->type!=NULL) { free(clean->type); clean->type=NULL; }
		    		if(clean->value!=NULL){ free(clean->value);clean->value=NULL;}
		    		free(clean);clean=NULL;
			}

		}//while
		head=NULL;
	}//if size
	if( goid != NULL) { free(goid); goid=NULL; }
}

int GO_find_id(struct GO *localhead, char **goid){
	struct GO *item=NULL;

	item = localhead;
	 while (item !=NULL){
		if( !strcmp(item->type,"GOid") ){
			strcpy(*goid, item->value);
			return 1;
		}
		item = item->next;
	}//while
	return 0; 
}
