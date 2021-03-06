%{
// $Id: parsedFunctionLexer.l,v 1.2 2002/12/12 10:22:09 oliver Exp $
#include <BALL/MATHS/parsedFunction.h>

using namespace BALL;

#include "parsedFunctionParser.h"
%}

%option noyywrap
%option nounput

LETTER   [a-zA-Z]
SIGN     "+"|"-"
DIGIT            [0-9]
INT      {SIGN}?{DIGIT}+
REAL     {SIGN}?({DIGIT}*)((\.({DIGIT}*))?)((E|e)({SIGN}?){DIGIT}{1,3})?
ID       {LETTER}(({LETTER}|{DIGIT})*)

%%
[ \t] /* Good bye white space...*/
"\n"    {return '\n';}
";"     {return ';';}
"\+"    {return  '+';}
"\-"    {return  '-';}
"\*"    {return  '*';}
"\/"    {return  '/';}
"\^"    {return  '^';}
"("     {return  '(';}
")"     {return  ')';}
"="     {return  '=';}
{INT}   {parsedFunctionlval.val=atoi(yytext); return NUM;}
{REAL}  {parsedFunctionlval.val=atof(yytext); return NUM;}
{ID}    {if (parsedFunctionFunctions->has(yytext)) 
				 {
					 double(*func)(double)=(*parsedFunctionFunctions)[yytext];
					 parsedFunctionlval.func = func;
					 return FNCT;
					}
					else
					{ if (parsedFunctionConstants->has(yytext))
					  {
							parsedFunctionlval.var = (*parsedFunctionConstants)[yytext];
							return VAR;
						}
						else
						{
							double *v = new double;
							(*parsedFunctionConstants)[yytext] = v;
							parsedFunctionlval.var = v;
							return VAR;
						}
					}
				}

%%
YY_BUFFER_STATE parsedFunction_buffer;
void parsedFunction_initBuffer(const char* buf)
{
	parsedFunction_buffer = parsedFunction_scan_string(buf);
}

void parsedFunction_delBuffer()
{
	parsedFunction_delete_buffer(parsedFunction_buffer);
}
