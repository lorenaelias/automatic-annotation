/* Compilador para separar dados pertinentes de termos do Gene Ontology */
%{
	#include <stdio.h>
	#include <string.h>
	#include <malloc.h>

	#define goidlist_maxlenght	10

	extern FILE *yyin;

	char **goidlist=NULL;

	char *lastgoid=NULL;

	int goidlist_index,idx;
	
	void yyerror(const char *str)
	{
		//fprintf(stderr,"error: %s\n",str);
		fprintf(stderr,".");
	}
	
	int yywrap()
	{
		return 1;
	}
	
	void goidlist_print (){
		for(idx=0; idx <= goidlist_index && idx < goidlist_maxlenght; idx++) 
			if (lastgoid != NULL && *(goidlist+idx) != NULL) {
				printf("INSERT INTO goidisa (goid, isa ) VALUES ('%s','%s');\n", lastgoid, *(goidlist+idx) );
				**(goidlist+idx)='\0';
			}
	}

	main( argc, argv )
	int argc;
	char **argv;
	{
	    ++argv, --argc;  /* skip over program name */

		goidlist = (char **)malloc( goidlist_maxlenght*sizeof(char *) );
		for( goidlist_index=0; goidlist_index<goidlist_maxlenght; goidlist_index++ ){
			*(goidlist+goidlist_index)=(char *)malloc( 20*sizeof(char) );
			**(goidlist+goidlist_index)='\0';
		}
		goidlist_index=-1;

		lastgoid = (char *)malloc( 20*sizeof(char) );
		*lastgoid = '\0';

		if ( argc > 0 )
			yyin = fopen( argv[0], "r" );
		else
			yyin = stdin;

		yyparse();
	}
%}

%union
{
	int integer;
	float real;
	char *string;
}

%token <string> ID
%token <string> NAME
%token <string> NAMESPACE
%token <string> DEF
%token <string> ECNUMBER
%token <string> ISA

%%

lines	: /* empty */
	| lines line
	;
	
line	: id name namespace def ecnumber isalist
			{ goidlist_print(); }
	|
	  id name namespace def isalist
	 { printf(",'');\n"); goidlist_print();}
	|
	  id name namespace isalist
	 { printf("'','');\n"); goidlist_print();}
	|	
	  error		
	;

isalist : /* empty list*/
	| isalist isa
	;

id 	: ID 
	  {	strcpy( lastgoid, strdup($1) );
			printf("INSERT INTO goterm (goid, \"name\", namespace, def, ecnumber) VALUES ('%s',", lastgoid); goidlist_index=-1;
		}
	;
	
name 	: NAME 
	  {printf("'%s',", $1); }
	;

namespace: NAMESPACE 
	  {printf("'%s',", $1); }
	;
	
def 	: DEF
	  {printf("'%s'", $1); }
	;

isa 	: ISA
	  {	goidlist_index++;
			if (goidlist_index < goidlist_maxlenght && goidlist_index >= 0 ) strncpy(*(goidlist+goidlist_index), strdup($1), 11); 
	  }
	;

ecnumber: ECNUMBER
	  {printf(",'%s');\n", $1);}
	 ;

%%

