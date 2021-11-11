/* Compilador para ler dados de um gene anotados em um arquivo EMBL */

%{
	#include <math.h>
	#include "gene.tab.h"
	extern YYSTYPE yylval;
	int size;

%}
FMARK	^"FT "
PSEUDO	"pseudo"
ID       [0-9a-zA-Z_:]*
CONTENT	["][ 0-9a-zA-Z_/:\-,%\n\.\(\);\*\+'\[\]`\|]+["]

%%

{FMARK}      

{PSEUDO}       	{	yylval.string=strdup(yytext);
			/*printf("id %s \n",yytext);*/
			return PSEUDO;
		}

{ID}        	{	yylval.string=strdup(yytext);
			/*printf("id %s \n",yytext);*/
			return IDENTIFIER;
		}

{CONTENT}      	{	
			yylval.string=strdup(yytext);
			/*printf("value %s \n",yytext);*/
			return VALUE;
		}

","	return COLUMN;

"%"	return PERCENT;

"-"	return MINUS;
		
"_"	return UND;		

";"	return SEMICOL;

"\""	return QUOTATION;

"*"	return STAR;

":"	return DOUBLEPOINT;

"."	return DOT;

"/"	return BAR;

"="	return EQUAL;

"("	return OPENPARENTHESIS;

")"	return CLOSEPARENTHESIS;

[ \t\r\n]+          /* eat up whitespace */

.           printf( "unknown %s\n", yytext );

%%

