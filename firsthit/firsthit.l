/* Compilador para separar dados pertinentes de termos do Gene Ontology */

%{
	#include <math.h>
	#include "firsthit.tab.h"
	extern YYSTYPE yylval;

	void substr(char *in, int from){
		int tam=strlen(in);
		int index=0;
		for (index=0; (index+from)<=tam; index++)
			in[index]=in[index+from];
		in[index+1]='\0';
		strcpy(yylval.string,in);
	}

%}

content           [^\n]*
content2          [^\]]*

%%

^"Query= "{content}     {	yylval.string=strdup(yytext);
				substr(yylval.string,7);
				return ID;
			}

^"         ("{content} 	{	yylval.string=strdup(yytext);
				substr(yylval.string,10);
				return TAMANHO_QUERY;
			}

"Length = "{content}   	{	yylval.string=strdup(yytext);
				substr(yylval.string,9);
				return TAMANHO_SUBJECT;
			}

^">"{content2}   	{	yylval.string=strdup(yytext);
				substr(yylval.string,1);
				return PRODUTO;
			}

"Expect = "{content}   	{	yylval.string=strdup(yytext);
				substr(yylval.string,9);
				return EVALUE;
			}
^" Identities = "{content}   	{	yylval.string=strdup(yytext);
					substr(yylval.string,14);
					return PERCENTUAL;
				}

[\t\r\n]+          /* eat up whitespace */

.	//printf("%s",yytext);

%%
