/* Compilador para separar dados pertinentes de termos do Gene Ontology */
%{
	#include <stdio.h>
	#include <string.h>
	#include <malloc.h>
	#define MAX_IDENTIFIER 300

	extern FILE *yyin;

	char splitted[9][MAX_IDENTIFIER];

	char similarto[6][MAX_IDENTIFIER];
	
	void splitdata(char *data, char *sep){
		int i;
		char *buffer=NULL;//(char *)malloc(MAX_IDENTIFIER);

		if(!data || !sep) return;
		buffer = strtok (data, sep);
		i=0;
		while (buffer!=NULL && i<10){
			memcpy(splitted[i],buffer,strlen(buffer)+1);
			buffer = strtok (NULL, sep);
			i=i+1;
		}	
	}
	
	void yyerror(const char *str)
	{
		//fprintf(stderr,"error: %s\n",str);
		fprintf(stderr,".");
	}
	
	int yywrap()
	{
		return 1;
	}
	
	main( argc, argv )
	int argc;
	char **argv;
	{
	    ++argv, --argc;  /* skip over program name */
	    
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
%token <string> TAMANHO_QUERY
%token <string> TAMANHO_SUBJECT
%token <string> PRODUTO
%token <string> EVALUE
%token <string> PERCENTUAL

%%

lines	: /* empty */
	| lines line
	;
	
line	: id tamanho_query produto tamanho_subject evalue percentual
	{printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\n",similarto[0], similarto[1], similarto[2], similarto[3], similarto[4], similarto[5], similarto[6]);}
	

	;
///similarity="Similar to <organismo>, <produto>, e-value: <valor>, <percentual>% id in <tamanho> aa"
id 	: ID 
	  { splitdata(strdup($1), " "); //printf("%s\t", splitted[0]);
		strcpy(similarto[0],splitted[0]);
	  }
	;
	
tamanho_query	: TAMANHO_QUERY
	  {splitdata(strdup($1), " "); //printf("%s\t", splitted[0]);
		strcpy(similarto[1],splitted[0]);
	  }
	;

produto: PRODUTO 
		//O produto possui formatos variados e sempre contem o nome ou sigla do produto proteico 
		//conjugado com o nome do organismo de origem. Foi necessario bolar regras diferentes para
		//cada situacao:
	  {splitdata(strdup($1), "|["); //Primeiro separa os vários pedaços de texto do produto
		if(!strcmp(splitted[0],"sp")){//Esse é o primeiro caso especial: 
			if(splitted[4][0]!='\0'){
				strcpy(similarto[2],splitted[4]);//produto proteico no quinto campo
				strcpy(similarto[3],splitted[5]);//organismo de origem no sexto campo
			}else{
				//Proteínas cujo prefixo é "sp" as vezes nem tem o nome do organismo de origem.
				strcpy(similarto[2],splitted[2]);
			}
		}else{//Se nao sao proteinas como prefixo sp entao fica mais facil e todas seguem
		//um padrao bem definido de disposicao do produto proteico e organismo de origem
		//Nesse caso o produto sempre esta na terceira posicao e organismo na quarta posicao.
			//Produto:
		    if(strchr(splitted[2],'\n')){ //Pode acontecer do produto estar quebrado em duas linhas ...
		   		splitdata(strdup(splitted[2]), "\n"); 
				strcpy(similarto[2],splitted[0]);//... e assim precisamos juntar esses dois pedacos ...
				strcat(similarto[2],splitted[1]);//... para recompor o produto proteico integral.
			}else{
				strcat(similarto[2],splitted[2]); //nao estao quebrado em duas linhas
	   		}
			//Organismo:
		    if(strchr(splitted[3],'\n')){//Pode acontecer do orgnanismo de origem estar quebrado em duas linhas ...
		   		splitdata(strdup(splitted[3]), "\n"); 
				strcpy(similarto[3],splitted[0]);//... e assim precisamos juntar esses dois pedacos ...
				strcat(similarto[3],splitted[1]);//... para recompor o produto proteico integral.
			}else{
				strcat(similarto[3],splitted[3]);//nao estao quebrado em duas linhas
	   		}
		}
	  }
//	| error
//	{ printf("No hits found\t-\t-\t-\n"); }
	;

tamanho_subject	: TAMANHO_SUBJECT
	  {splitdata(strdup($1), " "); //printf("%s\t", splitted[0]);
		strcpy(similarto[4],splitted[0]);
	  }
	;
	
evalue: EVALUE
	  {splitdata(strdup($1), ","); //printf("%s\t", splitted[0]); 
		strcpy(similarto[5],splitted[0]);
	  }
	 ;

percentual: PERCENTUAL
	  {splitdata(strdup($1), "()"); //printf("%s\n", splitted[1]); 
		strcpy(similarto[6],splitted[1]);
	  }
	 ;
%%

