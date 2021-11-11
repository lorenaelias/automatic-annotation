/* Compilador para separar dados pertinentes de termos do Gene Ontology */

%{
	#include <math.h>
	#include "go.tab.h"
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

%%

^"id: "{content}   	{	yylval.string=strdup(yytext);
				substr(yylval.string,4);
				return ID;
			}

^"name: "{content}   	{	yylval.string=strdup(yytext);
				substr(yylval.string,6);
				return NAME;
			}
^"namespace: "{content}   	{	yylval.string=strdup(yytext);
					substr(yylval.string,11);
					return NAMESPACE;
				}


^"def: "{content}   	{	yylval.string=strdup(yytext);
				substr(yylval.string,5);
				return DEF;
			}

^"xref: EC:"{content}   {	yylval.string=strdup(yytext);
				substr(yylval.string,6);
				return ECNUMBER;
			}

^"is_a: GO:"{content}   {	strncpy(yylval.string,strdup(yytext)+6,16);
				yylval.string[10]='\0';
				return ISA;
			}


[\t\r\n]+          /* eat up whitespace */

.

%%
