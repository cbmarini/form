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
int GetLinArgument(int *indexes, int *depth, mpf_t f_out, WORD *fun) {
	WORD *t, *tn, *tstop, *term, *arg, *argn, nargs, i;
	arg = fun + FUNHEAD;

	argn = arg; nargs = 0;
	while ( argn < fun + fun[1] ) {
		NEXTARG(argn);
		nargs++;
	}
	if ( nargs < 2 ) { 
		return(-1); 
	}
	/* Check that lin_ has two arguments*/
	if ( (*fun == LINFUNCTION) && (nargs != 2) ) { 
		return(-1); 
	}	

/* The first argument(s) should be small integer(s) */
	*depth = nargs - 1;
	for( i = 0; i < *depth; i++ ) {
		if (*arg != -SNUMBER) { return(-1); }
		indexes[i] = arg[1];
		NEXTARG(arg);
	}

/* The last argument can be a small number, fraction or float */
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
		return(0);
	}
	else if ( *t == FLOATFUN ) { /* float */
		UnpackFloat(f_out,t);
		if ( tn[-1] < 0 ) {/* change sign */
			mpf_neg(f_out,f_out);
		}
		return(0);
	}
	else { return(-1); }
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
	int numdigits = (AC.DefaultPrecision-AC.MaxWeight)*log10(2.0);
	GiNaC::Digits = numdigits;
/*
	We only consider -1 <= arg <= 1
	This can be changed later, but for now we don't burn ourselves 
	by picking a branch.
*/
	if ( (mpf_cmp_ui(arg,1L) > 0) || (mpf_cmp_si(arg,-1L) < 0) ) {
		return(-1);
	}
/*
	In principle this is well defined, so we can always implement this case. 
*/
	if ( weight < 1 ) {
		return(-1);
	}
/*
	Divergent case:	
*/
	if ( (weight == 1) && (mpf_cmp_ui(arg,1L) == 0) ) {
		return(-1);
	}
/*
	Convert to ginac numeric and evaluate
	Note: GiNaC evaluates lin_(1,-1) = 0????
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
  	#[ CalculateHpl :
*/
int CalculateHpl(mpf_t result, int *indexes, int depth, mpf_t arg) {
	int i;
	GiNaC::lst indexesg;
	char *str = nullptr;
	GiNaC::numeric argg;
	GiNaC::ex resultg;
	ostringstream oss;
/*
	Form measures precision in bits, GiNaC in decimal digits
*/
	int numdigits = (AC.DefaultPrecision-AC.MaxWeight)*log10(2.0);
	GiNaC::Digits = numdigits;
/*
	We only consider -1 <= arg <= 1
	This can be changed later, but for now we don't burn ourselves 
	by picking a branch.
*/
	if ( (mpf_cmp_ui(arg,1L) > 0) || (mpf_cmp_si(arg,-1L) < 0) ) {
		return(-1);
	}
/*
	Trailing zeroes lead to factors of hpl_(0;x) = ln_(x)
	So we only consider x > 0 in this case
*/
	if ( (indexes[depth-1] == 0) && (mpf_cmp_ui(arg,0L) <= 0) ) {
		return(-1);
	}
/*
	Leading ones lead to factors of hpl_(1;x) = ln_(1-x)
	So we only consider x < 1 in this case
*/
	if ( (indexes[0] == 1) && (mpf_cmp_ui(arg,1L) >= 0) ) {
		return(-1);
	}
/*
	Place the indexes in a ginac list
*/
	for( i=0; i < depth; i++ ) {
		indexesg.append(indexes[i]);
	}
/*
	Convert to ginac numeric and evaluate
*/
	gmp_asprintf(&str, "%.Ff", arg);
    argg = GiNaC::numeric(str);
	free(str);
	resultg = GiNaC::evalf(GiNaC::H(indexesg, argg));
/*
	Sometimes GiNaC returns a small complex number, even if the result is real.
*/
	resultg = GiNaC::real_part(resultg);
/*
	Convert back to mpf_t
*/
	oss << resultg;
	mpf_set_str(result, oss.str().c_str(), 10);
	return(0);
}
/*
  	#] CalculateHpl :
  	#[ EvaulatueLin : 
*/

int EvaluateLin(PHEAD WORD *term, WORD level, WORD par) {
	WORD *t, *tstop, *tt, *newterm, i, depth, first = 1;
	WORD *oldworkpointer = AT.WorkPointer, nsize, *indexes;
	int retval;

	DUMMYUSE(par);

	tstop = term + *term; tstop -= ABS(tstop[-1]);
	if ( AT.WorkPointer < term+*term ) AT.WorkPointer = term + *term;

/*
	Step 1: locate a LINFUNCTION
*/
	t = term+1;
	while ( t < tstop ) {
		if ( *t == LINFUNCTION || *t == HPLFUNCTION) {
			indexes = AT.WorkPointer;
			if (GetLinArgument(indexes,&depth,aux1,t) != 0) {
				MesPrint("Error: LINFUNCTION with illegal argument(s)");
				goto nextfun;
			}
/*
	Step 2: evaluate
*/
			if ( first ) {
				if ( *t == LINFUNCTION ) {
					if (CalculateLin(aux4,indexes[0],aux1) != 0) goto nextfun;
				}
				else if ( *t == HPLFUNCTION ) {
					if (CalculateHpl(aux4,indexes,depth,aux1) != 0) goto nextfun;
				}
				first = 0;
			}
			else {
				if ( *t == LINFUNCTION ) {
					if (CalculateLin(aux5,indexes[0],aux1) != 0) goto nextfun;
				}
				else if ( *t == HPLFUNCTION ) {
					/* HPL is not implemented yet */
					if (CalculateHpl(aux5,indexes,depth,aux1) != 0) goto nextfun;
				}
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
