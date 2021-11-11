/* Compilador para ler dados de um gene anotados em um arquivo EMBL */
%{
	#include <stdio.h>
	#include <string.h>
	#include <malloc.h>
	#include "GO.h"
	#define MAX_IDENTIFIER 300	
	extern FILE *yyin;
	char *feature_value,*buffer, *last_systematic_id;

	#define cds 1
	#define repetition 2
	#define domain 3
	#define trna 4	
	#define signal 5
	#define rrna 6
	#define gene 7		
	
	FILE *sql;

	struct type_feature{
		int type;
		int pos_begin;
		int pos_end;
		char orientation;
		char name[MAX_IDENTIFIER];
		char note[MAX_IDENTIFIER];
		char curation[MAX_IDENTIFIER];
		char predictedby[MAX_IDENTIFIER];
		char previous_systematic_id[MAX_IDENTIFIER];
		char product[MAX_IDENTIFIER];
		char domain_name[MAX_IDENTIFIER];
		char organism[MAX_IDENTIFIER];
		char EC_Number[MAX_IDENTIFIER];
		char systematic_id[MAX_IDENTIFIER];
		char protein_id[MAX_IDENTIFIER];
		char similarity[MAX_IDENTIFIER];
		char clevage[MAX_IDENTIFIER];
		char anticodon[MAX_IDENTIFIER];
		char blastp_file[MAX_IDENTIFIER];				
		int multipos[4][2];
		int pathogenicity;
		int pseudogene;
		int colour;
		char GeneID[MAX_IDENTIFIER];
		char GI[MAX_IDENTIFIER];
	};
	
	struct type_feature feature, feature_report;
	
	void featureinit(int previos_pos_begin, int previous_pos_end){
		int i;

		if ( 	feature.pos_begin != previos_pos_begin || 
					feature.pos_end   != previous_pos_end     ){
			feature.pos_begin=0;
			feature.pos_end=0;
			feature.previous_systematic_id[0]='\0';
			feature.systematic_id[0]='\0';
		}
		feature.type= 0;	
		feature.orientation='\0';
		feature.name[0]='\0';
		feature.note[0]='\0';
		feature.curation[0]='\0';
		feature.predictedby[0]='\0';
		feature.product[0]='\0';
		feature.domain_name[0]='\0';
		feature.organism[0]='\0';
		feature.EC_Number[0]='\0';
		feature.protein_id[0]='\0';
		feature.similarity[0]='\0';
		feature.clevage[0]='\0';
		feature.anticodon[0]='\0';
		feature.blastp_file[0]='\0';
		for(i=0; i<=4; feature.multipos[i][0]=0,feature.multipos[i][1]=0,i++);
		feature.pathogenicity=0;
		feature.pseudogene=0;
		feature.colour=3;
		feature.GI[0]='\0';
		feature.GeneID[0]='\0';
	}

	char *splitted[9];
	
	void splittedinit(){
		int i;
		for(i=0; i<10; i++){
			if(!splitted[i]) 
				splitted[i]=(char *)malloc(MAX_IDENTIFIER);
			splitted[i][0]='\0';
		}
	}
	
	void splitdata(char *data, char *sep){
		int i;
		char *buffer=NULL;//(char *)malloc(MAX_IDENTIFIER);
		
		if(!data || !sep) return;
		splittedinit();
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
		/*fprintf(stderr,"error: %s\n",str);*/
		fprintf(stderr,".");
		
	}
	
	int yywrap()
	{
		return 1;
	}
	
	void remove_embl_longtextformat(char *in, char *out){
		int	 pos;
		char *buffer=(char *)malloc(MAX_IDENTIFIER-10);
		out[0]='\0';
		/*Somente remove os caracteres indesejáveis se eles existirem*/
		if(strchr(in,'\n')!=NULL){
		
			/*	Primeiro remove ocorrencias indicando um qualifier na linha do arquivo embl.
				Todas as linhas de qualifier sempre começam com o mnemonico FT seguido de 
				19 espaços em branco, sempre. Por isso o texto abaixo tem esse formato.		*/
			do{				
				/*procura onde esta o primeira posicao de ocorrencia do FT*/
				buffer = strstr(in, "FT                   ");
				/*se encontrou uma posicao de ocorrencia entao prossiga */
				if(buffer!=NULL){
					/* calcula a posicao inicial do corte da sequencia FT*/
					pos=strlen(in)-strlen(buffer);	/*posição alvo*/
					/*	move as 21 posicoes que ficam na frente da sequencia FT 
						para cima das 21 posicoes que contem a sequencia FT, ou
						seja, a sequencia FT esta sendo sobrescrita com conteudo 
						valido
					*/
					memmove(in+pos,in+pos+21,strlen(in)-21);
				}
				/*quando nao achar mais posicoes FT entao pare */
			}while(buffer!=NULL);
			
		
			/*Depois substitui as quebras de linhas por espaços em branco*/
			buffer = strtok (in,"\n");
			while (buffer != NULL){
				strcat(out,strdup(buffer));
				buffer = strtok (NULL, "\n");
				if(buffer!=NULL) strcat(out," ");
			}	
			
		}else{ 
			strcpy(out,in);
		}
	}	

	void remove_character(char *in, char *out, char *undesirable){
		char *buffer=NULL;
		out[0]='\0';
		/*Somente remove o caracter indesejável se ele existir*/
		if(strchr(in, *undesirable)!=NULL){
			buffer = strtok (in,undesirable);
			while (buffer != NULL){
				if(strchr(in, undesirable)==NULL){
					strcat(out,strdup(buffer));
				}
				buffer = strtok (NULL, undesirable);
			}
		}else{ 
			strcpy(out,in);
		}
	}	

	
	void ImprimirUmRegistro(struct type_feature feature_report){
		int index;
	  /* Imprime o relatório da ultima feature carregada */
	    if(feature_report.type!=0)
		  switch(feature_report.type){
		  case gene:
		  case cds:
			if(feature_report.type == gene && feature_report.pseudogene == 0 ){
			/* Quando existe uma feature gene e uma feature CDS com o mesmo locus_tag devemos armazenar somente a entidade CDS no banco */
				return;
			}
		  	sql=fopen("./insert.gene.sql","a");
		  	if(sql!=NULL){
			  //fprintf(sql,"update GENE set blastp_file='%s' where systematic_id='%s';\n", feature_report.blastp_file, feature_report.systematic_id);
		  	  fprintf(sql,"insert into GENE (systematic_id, pos_begin, pos_end, orientation, name, protein_id,");
		  	  fprintf(sql,"curation, predictedby, previous_systematic_id, product, organism, EC_Number, ");
		  	  fprintf(sql,"similarity, blastp_file, pseudogene, pathogenicity, colour, type) ");		  	  
//		  	  fprintf(sql,"similarity, blastp_file, pseudogene, pathogenicity, colour, type, GI, GeneID) ");		  	  
	  	  fprintf(sql," values ('%s', %d, %d, '%c', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', %d, %d, %d, %d);\n",
//	  	  fprintf(sql," values ('%s', %d, %d, '%c', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', %d, %d, %d, %d, '%s', '%s');\n",
		  	  feature_report.systematic_id, feature_report.pos_begin, feature_report.pos_end,
		  	  feature_report.orientation, feature_report.name, feature_report.protein_id, feature_report.curation,
		  	  feature_report.predictedby, feature_report.previous_systematic_id, feature_report.product,
		  	  feature_report.organism, feature_report.EC_Number, feature_report.similarity,
		  	  feature_report.blastp_file, feature_report.pseudogene, feature_report.pathogenicity, feature_report.colour, 
			  feature_report.type); //, feature_report.GI, feature_report.GeneID);
		  	  fclose(sql);
		  	}
		  	sql=fopen("./insert.GO.sql","a");
		  	if(sql!=NULL){
		  		GO_report(sql, feature_report.systematic_id);
		  		fclose(sql);
		  	}
			if (feature_report.multipos[0][0] > 0){
			   sql=fopen("./insert.multipos.sql","a");
			   if(sql!=NULL){
			   	for(index=0; feature_report.multipos[index][0]; index ++){
				  fprintf(sql,"insert into MULTIPOS (GENE_systematic_id, pos_begin, pos_end) values ('%s', %d, %d);\n",
				  feature_report.systematic_id, feature_report.multipos[index][0], feature_report.multipos[index][1]);
				}
			  }
			  fclose(sql);
			}
		  	break;
		  case repetition:
		  	sql=fopen("./insert.repeat.sql","a");
		  	if(sql!=NULL){
		  	  fprintf(sql,"insert into  DNAREPEAT (identifiedby, pos_begin, pos_end,orientation, note) ");
		  	  fprintf(sql," values ('%s', %d, %d, '%c', '%s');\n",
		  	  feature_report.predictedby, feature_report.pos_begin, feature_report.pos_end,
		  	  feature_report.orientation, feature_report.note);
		  	  fclose(sql);
		  	}
		  	break;
		  case domain:
		  	sql=fopen("./insert.domain.sql","a");
		  	if(sql!=NULL){
		  	  fprintf(sql,"insert into  DOMAIN (GENE_systematic_id, label, pos_begin, pos_end,orientation, previous_systematic_id, name, note) ");
		  	  fprintf(sql," values ('%s', '%s', %d, %d, '%c', '%s', '%s', '%s');\n",
		  	  last_systematic_id, feature_report.predictedby, feature_report.pos_begin, feature_report.pos_end,
		  	  feature_report.orientation, feature_report.previous_systematic_id, feature_report.domain_name,
		  	  feature_report.note);
		  	  fclose(sql);
		  	}
		  	break;
		  case signal:
		  	sql=fopen("./insert.signal.sql","a");
		  	
		  	/* Matriz splitted recebeu dados separados por virgula */
		  	splitdata(feature_report.clevage, ","); 
		  	
		  	if(sql!=NULL){
		  	  fprintf(sql,"insert into  SIGNAL (GENE_systematic_id, predictedby, pos_begin, pos_end,orientation, previous_systematic_id, pos_clevage, signaltext, note) ");
		  	  fprintf(sql," values ('%s', '%s', %d, %d, '%c', '%s', %s, '%s', '%s');\n",
		  	  last_systematic_id, feature_report.predictedby, feature_report.pos_begin, feature_report.pos_end,
		  	  feature_report.orientation, feature_report.previous_systematic_id, splitted[0],
		  	  feature_report.product, feature_report.note);
		  	  fclose(sql);
		  	}
		  	break;
		  case trna:
		  	sql=fopen("./insert.trna.sql","a");
		  	
		  	if(sql!=NULL){
		  	  fprintf(sql,"insert into  tRNA (locus_tag, pos_begin, pos_end, orientation, product, note, anticodon, aa, organism) ");
		  	  fprintf(sql," values ('%s', %d, %d, '%c', '%s', '%s', '%s', '%s', '%s');\n",
		  	  feature_report.systematic_id, feature_report.pos_begin, feature_report.pos_end, feature_report.orientation, 
			  feature_report.product, feature_report.note, feature_report.anticodon, feature_report.name, feature_report.organism);
		  	  fclose(sql);
		  	}
		  	break;
		  case rrna:
		  	sql=fopen("./insert.rrna.sql","a");
		  	
		  	if(sql!=NULL){
		  	  fprintf(sql,"insert into  rRNA (locus_tag, pos_begin, pos_end, orientation, product, organism) ");
		  	  fprintf(sql," values ('%s', %d, %d, '%c', '%s', '%s');\n",
		  	  feature_report.systematic_id, feature_report.pos_begin, feature_report.pos_end,
		  	  feature_report.orientation, feature_report.product, feature_report.organism);
		  	  fclose(sql);
		  	}
		  	break;
		  	
		  }
	}/*fim GerarRelatorio */
	
	main( argc, argv )
	int argc;
	char **argv;
	{
	    ++argv, --argc;  /* skip over program name */
	    
	    /* Initialize buffer to format data extracted by parser*/
	    feature_value = (char *)malloc(MAX_IDENTIFIER);
	    feature_value[0]='\0';
	    buffer = (char *)malloc(MAX_IDENTIFIER);
	    buffer[0]='\0';
	    last_systematic_id= (char *)malloc(MAX_IDENTIFIER);
	    last_systematic_id[0]='\0';
	    
	    
	    if ( argc > 0 )
		    yyin = fopen( argv[0], "r" );
	    else
		    yyin = stdin;

		
	    featureinit(-1,-1);
	    feature_report = feature;

	    yyparse();
	    /* 	Imprime o ultimo registro porque nao ha mais novos 
	    	registros para iniciar a variavel 'feature_report' */
	    ImprimirUmRegistro(feature); 
	    printf("\n\nReference:\nSantos, A R; Silva A; Miyoshi, Azevedo, V. (2011) CpDB: A relational database schema and tools for bacterial genomes annotation and pos-genome research (unpublished).\nPlease, contact me before use or cite it for eventual updates: 'asantos.icb.ufmg.br@gmail.com'.\nThank you, Anderson Santos.\n");
	}
%}

%union
{
	int integer;
	float real;
	char *string;
}

%token <integer> POSITION 
%token <real>	QUALITY
%token <string> IDENTIFIER
%token <string> VALUE
%token <string> PSEUDO
%token COLUMN
%token PERCENT
%token MINUS
%token UND
%token SEMICOL
%token QUOTATION
%token STAR
%token DOUBLEPOINT
%token DOT
%token BAR
%token EQUAL
%token OPENPARENTHESIS
%token CLOSEPARENTHESIS
%token FMARK

%%

lines	: /* empty */
	| lines line
	;
	
line	:	IDENTIFIER coords
		{ printf("Tipo de feature '%s'\n", $1);

		  if(!strcmp($1,"CDS")){
			feature.type= cds;
		  }else if(!strcmp($1,"gene")){
			feature.type= gene;
		  }else if(!strcmp($1,"repeat_region")){
			feature.type= repetition;
		  }else if(!strcmp($1,"misc_feature")){
			feature.type= domain;
			/* registros de peptideo sinal sao diferenciados de dominios pela 
			   existencia de um qualificador denominado 'coord', mas
			   inicialmente sao denominados pertencentes ao tipo dominio */
		  }else if(!strcmp($1,"tRNA")){
			feature.type= trna;
		  }else if(!strcmp($1,"rRNA")){
			feature.type= rrna;
		  }
		  /* A cada nova entrada encontrada faz a impressao do registro anteriormente lido */	
		  ImprimirUmRegistro(feature_report);
		  
		  strcpy(feature.organism, "Corynebacterium pseudoturbeculosis");
		  		  
		}
	|	IDENTIFIER EQUAL VALUE
		{	
			//O parser não suporta o campo translation. Assim ele é evitado.
			if(strcmp($1,"translation")!=0 && strcmp($1,"note")!=0 ){
				/* Antes de fazer uso do valor da feature remove as aspas
				 do inicio e do fim da string que possui o valor desejado */
				strcpy(buffer,strdup($3));
				remove_character(strdup($3), buffer, "\"");
				remove_embl_longtextformat(strdup(buffer), feature_value);
				printf("Feature '%s' com valor '%s'\n", $1, feature_value);
			}
			if(!strcmp($1,"curation")){
				strcpy(feature.curation, feature_value);
			}else if(!strcmp($1,"previous_systematic_id") || !strcmp($1,"id")){
				strcpy(feature.previous_systematic_id, feature_value);
			}else if(!strcmp($1,"systematic_id") || !strcmp($1,"locus_tag")){
				strcpy(feature.systematic_id, feature_value);
				if(last_systematic_id!=NULL){
					last_systematic_id[0]='\0';
					strcpy(last_systematic_id, strdup(feature_value));
				}
			}else if(!strcmp($1,"product")){
				strcpy(feature.product, feature_value);
			}else if(!strcmp($1,"protein_id")){
				strcpy(feature.protein_id, feature_value);
			}else if(!strcmp($1,"domain")){
				strcpy(feature.domain_name, feature_value);
			}else if(!strcmp($1,"EC_Number")){
				strcpy(feature.EC_Number, feature_value);
			}else if(!strcmp($1,"gene")){
				strcpy(feature.name, feature_value);
			}else if(!strcmp($1,"similarity")){
				strcpy(feature.similarity, feature_value);
			}else if(!strcmp($1,"repeat_region") || !strcmp($1,"label")){
				strcpy(feature.predictedby, feature_value);
			}else if(!strcmp($1,"sig_peptide")){
				feature.type = signal;
				strcpy(feature.product, feature_value);
			}else if(!strcmp($1,"blastp_file")){
				strcpy(feature.blastp_file, feature_value);
			}else if(!strcmp($1,"coord")){
				feature.type = signal;
				strcpy(feature.clevage, feature_value);
			}else if(!strcmp($1,"anticodon")){
				strcpy(feature.anticodon, feature_value);
			}else if(!strcmp($1,"GO_component")){
				//GO_add( strdup("component"), strdup(feature_value) );
			}else if(!strcmp($1,"GO_function")){
				//GO_add( strdup("function"), strdup(feature_value) );
			}else if(!strcmp($1,"GO_process")){
				//GO_add( strdup("process"), strdup(feature_value) );
			}else if(!strcmp($1,"pathogenicity")){
				feature.pathogenicity = atoi(strdup(feature_value));
			}else if(!strcmp($1,"note")){
				if(strstr(feature_value,"Glimmer")){
					strcpy(feature.predictedby, "glimmer");	
				}else if(strstr(feature_value,"manually")){
					strcpy(feature.predictedby, "human");	
				}else	
					strcpy(feature.note, feature_value);	 
			}else if(!strcmp($1,"db_xref")){
				if(strstr(feature_value,"GI:")){
					strncpy(feature.GI, feature_value+3, strlen(feature_value));	
				}else if(strstr(feature_value,"GeneID:")){
					strncpy(feature.GeneID, feature_value+7, strlen(feature_value));	
				}else	
					strcpy(feature.note, feature_value);	 
			}
		}
		
	|	IDENTIFIER EQUAL IDENTIFIER
		{	
			//O parser não suporta o campo translation. Assim ele é evitado.
			if(strcmp($1,"translation")!=0 && strcmp($1,"note")!=0 ){
				remove_embl_longtextformat(strdup($3), feature_value);
			}
			printf("Feature '%s' com valor '%s'\n", $1, feature_value);
			if(!strcmp($1,"label")){
				strcpy(feature.predictedby, feature_value);
			}else if(!strcmp($1,"GO:aspect") || !strcmp($1,"aspect")){
				GO_add( strdup("begin"), strdup("begin") );
				GO_add( strdup("aspect"), strdup(feature_value) );
			}else if(!strcmp($1,"term")){
				GO_add( strdup("term"), strdup(feature_value) );
			}else if(!strcmp($1,"GOid")){
				GO_add( strdup("GOid"), strdup(feature_value) );
			}else if(!strcmp($1,"evidence")){
				GO_add( strdup("evidence"), strdup(feature_value) );
			}else if(!strcmp($1,"colour")){
				feature.colour = atoi(strdup(feature_value));			
			}else if(!strcmp($1,"pathogenicity")){
				feature.pathogenicity = atoi(strdup(feature_value));
			}
		}

	|	IDENTIFIER EQUAL OPENPARENTHESIS 
					IDENTIFIER DOUBLEPOINT IDENTIFIER DOT DOT IDENTIFIER COLUMN 
					IDENTIFIER DOUBLEPOINT IDENTIFIER 
				 CLOSEPARENTHESIS
		{	
			if(!strcmp($1,"anticodon")){
				sprintf(feature.anticodon,"pos(%s,%s), aa: %s", $6, $9, $13);
				strcpy(feature.name,$13);
				printf("Feature '%s' com valor '%s'\n", $1, feature.anticodon);		
			}
		}

	|	PSEUDO
		{
			printf("Feature '%s' com valor '%d'\n", $1, feature.pseudogene);
			feature.pseudogene = 1;
		}

	|	error
	;

coords: 	IDENTIFIER DOT DOT IDENTIFIER
		{
			/* 
				Esse eh um ponto critico do programa. Quando entramos na regra coords
				estamos comecando a capturar dados de uma nova feature e caso existam dados
				carregados na estrutura 'gene' sobre a ultima feature lida precisamos descarregar
				esses dados antes que a estrutura de armazenamento temporario 'gene' seja limpa
				para uma nova rodada de carregamento de dados. Assim eh preciso salvar os dados lidos
				anteriormente para impressao no momento que a regra coords deixar de ser avaliada.
				Nesse caso a impressao dos dados carregados anteriormente acontecera no fim da regra
				'line', que foi a chamadora da regra 'coords'.
			*/
			 printf("\ndireta  ini '%s', fim '%s'\n", $1, $4);
			 feature_report = feature; featureinit(atoi($1), atoi($4));
			 feature.orientation='+'; feature.pos_begin=atoi($1); feature.pos_end=atoi($4);
		}
	|	IDENTIFIER OPENPARENTHESIS IDENTIFIER DOT DOT IDENTIFIER COLUMN 
					   IDENTIFIER DOT DOT IDENTIFIER 	CLOSEPARENTHESIS
		{
			 printf("\ndireta pseudo2  ini '%s', fim '%s'\n", $3, $11);
			 feature_report = feature; featureinit(atoi($3), atoi($11));
			 feature.orientation='+'; feature.pos_begin=atoi($3); feature.pos_end=atoi($11);
			 feature.multipos[0][0]=atoi($3);feature.multipos[0][1]=atoi($6);
			 feature.multipos[1][0]=atoi($8);feature.multipos[1][1]=atoi($11);
		}		
	|	IDENTIFIER OPENPARENTHESIS IDENTIFIER DOT DOT IDENTIFIER COLUMN 
					   IDENTIFIER DOT DOT IDENTIFIER COLUMN 
					   IDENTIFIER DOT DOT IDENTIFIER 	CLOSEPARENTHESIS
		{
			 printf("\ndireta pseudo3  ini '%s', fim '%s'\n", $3, $16);
			 feature_report = feature; featureinit(atoi($3), atoi($16));
			 feature.orientation='+'; feature.pos_begin=atoi($3); feature.pos_end=atoi($16);
			 feature.multipos[0][0]=atoi($3 );feature.multipos[0][1]=atoi($6);
			 feature.multipos[1][0]=atoi($8 );feature.multipos[1][1]=atoi($11);
			 feature.multipos[2][0]=atoi($13);feature.multipos[2][1]=atoi($16);
		}		
	|	IDENTIFIER OPENPARENTHESIS IDENTIFIER DOT DOT IDENTIFIER COLUMN 
					   IDENTIFIER DOT DOT IDENTIFIER COLUMN 
					   IDENTIFIER DOT DOT IDENTIFIER COLUMN 
					   IDENTIFIER DOT DOT IDENTIFIER 	CLOSEPARENTHESIS
		{
			 printf("\ndireta pseudo4  ini '%s', fim '%s'\n", $3, $21);
			 feature_report = feature; featureinit(atoi($3), atoi($21));
			 feature.orientation='+'; feature.pos_begin=atoi($3); feature.pos_end=atoi($21);
			 feature.multipos[0][0]=atoi($3 );feature.multipos[0][1]=atoi($6);
			 feature.multipos[1][0]=atoi($8 );feature.multipos[1][1]=atoi($11);
			 feature.multipos[2][0]=atoi($13);feature.multipos[2][1]=atoi($16);
			 feature.multipos[3][0]=atoi($18);feature.multipos[3][1]=atoi($21);
		}		
	|	IDENTIFIER OPENPARENTHESIS IDENTIFIER DOT DOT IDENTIFIER CLOSEPARENTHESIS
		{
			 printf("\nreversa ini '%s', fim '%s'\n", $3, $6);
			 feature_report = feature; featureinit(atoi($3), atoi($6));
			 feature.orientation='-'; feature.pos_begin=atoi($3); feature.pos_end=atoi($6);
		}
	|	IDENTIFIER	OPENPARENTHESIS 
			IDENTIFIER	OPENPARENTHESIS 
				IDENTIFIER DOT DOT IDENTIFIER COLUMN
				IDENTIFIER DOT DOT IDENTIFIER
					CLOSEPARENTHESIS
				CLOSEPARENTHESIS							
		{
			 printf("\nreversa pseudo2 ini '%s', fim '%s'\n", $5, $13);
			 feature_report = feature; featureinit(atoi($5), atoi($13));
			 feature.orientation='-'; feature.pos_begin=atoi($5); feature.pos_end=atoi($13);
			 feature.multipos[0][0]=atoi($5 );feature.multipos[0][1]=atoi($8);
			 feature.multipos[1][0]=atoi($10);feature.multipos[1][1]=atoi($13);
		}
	|	IDENTIFIER	OPENPARENTHESIS 
			IDENTIFIER	OPENPARENTHESIS 
				IDENTIFIER DOT DOT IDENTIFIER COLUMN
				IDENTIFIER DOT DOT IDENTIFIER COLUMN				
				IDENTIFIER DOT DOT IDENTIFIER
					CLOSEPARENTHESIS
				CLOSEPARENTHESIS							
		{
			 printf("\nreversa pseudo3 ini '%s', fim '%s'\n", $5, $18);
			 feature_report = feature; featureinit(atoi($5), atoi($18));
			 feature.orientation='-'; feature.pos_begin=atoi($5); feature.pos_end=atoi($18);
			 feature.multipos[0][0]=atoi($5 );feature.multipos[0][1]=atoi($8);
			 feature.multipos[1][0]=atoi($10);feature.multipos[1][1]=atoi($13);
			 feature.multipos[2][0]=atoi($15);feature.multipos[2][1]=atoi($18);
		}
	|	IDENTIFIER	OPENPARENTHESIS 
			IDENTIFIER	OPENPARENTHESIS 
				IDENTIFIER DOT DOT IDENTIFIER COLUMN
				IDENTIFIER DOT DOT IDENTIFIER COLUMN				
				IDENTIFIER DOT DOT IDENTIFIER COLUMN				
				IDENTIFIER DOT DOT IDENTIFIER
					CLOSEPARENTHESIS
				CLOSEPARENTHESIS							
		{
			 printf("\nreversa pseudo4 ini '%s', fim '%s'\n", $5, $23);
			 feature_report = feature; featureinit(atoi($5), atoi($23));
			 feature.orientation='-'; feature.pos_begin=atoi($5); feature.pos_end=atoi($23);
			 feature.multipos[0][0]=atoi($5 );feature.multipos[0][1]=atoi($8);
			 feature.multipos[1][0]=atoi($10);feature.multipos[1][1]=atoi($13);
			 feature.multipos[2][0]=atoi($15);feature.multipos[2][1]=atoi($18);
			 feature.multipos[3][0]=atoi($20);feature.multipos[3][1]=atoi($23);
		}
	;

%%
