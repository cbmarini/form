/* 
  	#[ Includes :
*/
#include <iostream>
#include <sstream>
#include <ginac/ginac.h>
#include "gmp.h"
using namespace std;

extern "C" {
/* First include the header files of ginac to avoid macro collision */
	#include "form3.h"
	int PackFloat(WORD *,mpf_t);
	int UnpackFloat(mpf_t,WORD *);
	void RatToFloat(mpf_t, UWORD *, int);
	void FormtoZ(mpz_t,UWORD *,WORD);
}
/*
  	#] Includes :
  	#[ GetLinArgument :
*/
int GetLinArgument(int *weight, mpf_t f_out, WORD *fun) {
	WORD *t, *tn, *tstop, *term, *arg, *argn, ncoef;
	arg = fun + FUNHEAD;
/* Check that lin_ has two arguments*/
	argn = arg; ncoef = 0;
	while ( argn < fun + fun[1] ) {
		NEXTARG(argn);
		ncoef++;
	}
	if ( ncoef != 2) { 
		return(-1); 
	}	

/* The first argument should be a small integer */
	if (*arg != -SNUMBER) { 
		return(-1); 
	}
	*weight = arg[1];

/* The second argument can be a small number, fraction or float */
	NEXTARG(arg);
	if (*arg == -SNUMBER) { /* small number */
		mpf_set_si(f_out,arg[1]);
		return(0);
	}

	term = arg + ARGHEAD;
	tn = term + *term;
	tstop = tn - ABS(tn[-1]);
	t = term + 1;
	if ( t == tstop) { /* fraction */
		RatToFloat(f_out,(UWORD *)t,tn[-1]);
	}
	else if ( *t == FLOATFUN ) { /* float */
		UnpackFloat(f_out,t);
		if ( tn[-1] < 0 ) {/* change sign */
			mpf_neg(f_out,f_out);
		}
	}
	else { return(-1); }

	return(0);
}
/*
  	#] GetLinArgument :
  	#[ CalculateLin :
*/
int CalculateLin(mpf_t result, int weight, mpf_t arg) {
	char *str = nullptr;
	GiNaC::numeric argg;
	GiNaC::ex resultg;
	ostringstream oss;
/*
	Form measures precision in bits, GiNaC in decimal digits
*/
	int numdigits = (AC.DefaultPrecision-AC.MaxWeight-1)*log10(2.0);
	GiNaC::Digits = numdigits;
/*
	We only consider -1 <= arg <= 1:
*/
	if ( (mpf_cmp_ui(arg,1L) > 0) || (mpf_cmp_si(arg,-1L) < 0) ) {
		return(1);
	}
/*
	Convert to ginac numeric and evaluate
*/
	gmp_asprintf(&str, "%.Ff", arg);
    argg = GiNaC::numeric(str);
	free(str);
	resultg = GiNaC::evalf(GiNaC::Li(weight, argg));
/*
	Convert back to mpf_t
*/
	oss << resultg;
	mpf_set_str(result, oss.str().c_str(), 10);
	return(0);
}
/*
  	#] CalculateLin :
  	#[ EvaulatueLin : 
*/

int EvaluateLin(PHEAD WORD *term, WORD level, WORD par) {
	WORD *t, *tstop, *tt, *newterm, i, weight, first = 1;
	WORD *oldworkpointer = AT.WorkPointer, nsize;
	int retval;

	DUMMYUSE(par);

	tstop = term + *term; tstop -= ABS(tstop[-1]);
	if ( AT.WorkPointer < term+*term ) AT.WorkPointer = term + *term;

/*
	Step 1: locate a LINFUNCTION
*/
	t = term+1;
	while ( t < tstop ) {
		if ( *t == LINFUNCTION ) {
			if (GetLinArgument(&weight,aux1,t) != 0) {
				MesPrint("Error: LINFUNCTION with illegal argument(s)");
				goto nextfun;
			}
/*
	Step 2: evaluate
*/
			if ( first ) {
				if (CalculateLin(aux4,weight,aux1) != 0) goto nextfun;
				first = 0;
			}
			else {
				if (CalculateLin(aux5,weight,aux1) != 0) goto nextfun;
				mpf_mul(aux4,aux4,aux5);
			}
			*t = 0;
		}
nextfun:
		t += t[1];
	}
	if ( first == 1 ) return(Generator(BHEAD term,level));

/*
	Step 3:
	Now the regular coefficient, if it is not 1/1.
	We have two cases: size +- 3, or bigger.
*/
	nsize = term[*term-1];
	if ( nsize < 0 ) {
		mpf_neg(aux4,aux4);
		nsize = -nsize;
	}
	if ( nsize == 3 ) {
		if ( tstop[0] != 1 ) {
			mpf_mul_ui(aux4,aux4,(ULONG)((UWORD)tstop[0]));
		}
		if ( tstop[1] != 1 ) {
			mpf_div_ui(aux4,aux4,(ULONG)((UWORD)tstop[1]));
		}
	}
	else {
		RatToFloat(aux5,(UWORD *)tstop,nsize);
		mpf_mul(aux4,aux4,aux5);
	}
/*
	Now we have to locate possible other float_ functions.
*/
	t = term+1;
	while ( t < tstop ) {
		if ( *t == FLOATFUN ) {
			UnpackFloat(aux5,t);
			mpf_mul(aux4,aux4,aux5);
		}
		t += t[1];
	}

/*
	Now we should compose the new term in the WorkSpace.
*/
	t = term+1;
	newterm = AT.WorkPointer;
	tt = newterm+1;
	while ( t < tstop ) {
		if ( *t == 0 || *t == FLOATFUN ) t += t[1];
		else {
			i = t[1]; NCOPY(tt,t,i);
		}
	}
	PackFloat(tt,aux4);
	tt += tt[1];
	*tt++ = 1; *tt++ = 1; *tt++ = 3;
	*newterm = tt-newterm;
	AT.WorkPointer = tt;
	retval = Generator(BHEAD newterm,level);
	AT.WorkPointer = oldworkpointer;
	return(retval);
}
/*
  	#] EvaluateLin :
*/
