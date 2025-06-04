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
  	#[ ParseArgumentToFloat :
*/
int ParseArgumentToFloat(mpf_t f_out, WORD *arg) {
    WORD *term, *tn, *tstop, *t;

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

    return -1;  // unsupported format
}

/*
  	#] ParseArgumentToFloat
  	#[ GetMplArgument :
*/
int GetMplArgument(int *indexes, int *depth, mpf_t *f_out, WORD *fun) {
	WORD *arg, *argnext, *argstop, *term, *tn, *tstop, *t, i;
/* 
	Get indexes from the first argument
*/
	arg = fun + FUNHEAD;
	// Check that the argument is not a fast argument
	if ( *arg < 0 ) { 
		return(-1); 
	}
	argnext = arg+*arg;
	argstop = argnext - ABS(argnext[-1]);
	// Check that the argument is a lst_ function that is properly normalized
	if ( (argnext[-3] != 1) || (argnext[-2] != 1) || (argnext[-1] != 3) ) { 
		return(-1); 
	}
	arg = arg + ARGHEAD; arg = arg + 1;
	if ( *arg != LSTFUNCTION || (arg+arg[1]) != argstop ) { 
		return(-1); 
	}
	// We can now read the indexes
	arg = arg + FUNHEAD;
	*depth = 0;
	while ( arg < argstop ) {
		if ( *arg != -SNUMBER ) return(-1); 
		indexes[*depth] = arg[1];
		(*depth)++;
		NEXTARG(arg);
	}
/* 
	Get the arguments from the second argument
*/
	arg = fun + FUNHEAD;
	NEXTARG(arg);
	// Check that the argument is not a fast argument
	if ( *arg < 0 ) { 
		return(-1); 
	}
	argnext = arg+*arg;
	argstop = argnext - ABS(argnext[-1]);
	// Check that the argument is a lst_ function that is properly normalized
	if ( (argnext[-3] != 1) || (argnext[-2] != 1) || (argnext[-1] != 3) ) { 
		return(-1); 
	}
	arg = arg + ARGHEAD; arg = arg + 1;
	if ( *arg != LSTFUNCTION || (arg+arg[1]) != argstop ) { 
		return(-1); 
	}
	// We can now read the arguments
	arg = arg + FUNHEAD;
	i = 0;
	while( arg < argstop ) {
		if ( ParseArgumentToFloat(f_out[i], arg) != 0 ) {
			return(-1); 
		}
		i++; NEXTARG(arg);
	}
	// Check that we have the right number of arguments
	if ( i != *depth ) { 
		return(-1); 
	}
	return(0);
}
/*
  	#] GetMplArgument :
  	#[ GetLinArgument :
*/
int GetLinArgument(int *indexes, int *depth, mpf_t *f_out, WORD *fun) {
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
	if ( (*fun == LINFUNCTION || *fun == MPLFUNCTION) && (nargs != 2) ) { 
		return(-1); 
	}	
	if ( *fun == MPLFUNCTION ) {
		return(GetMplArgument(indexes, depth, f_out, fun));
	}
/* The first argument(s) should be small integer(s) */
	*depth = nargs - 1;
	for( i = 0; i < *depth; i++ ) {
		if (*arg != -SNUMBER) { return(-1); }
		indexes[i] = arg[1];
		NEXTARG(arg);
	}

/* The last argument can be a small number, fraction or float */
	if (ParseArgumentToFloat(f_out[0], arg) != 0) {
		return(-1);
	}
	return(0);
}
/*
  	#] GetLinArgument :
  	#[ CalculateLin :
*/
int CalculateLin(mpf_t resultRe, mpf_t resultIm, int weight, mpf_t arg) {
	char *str = nullptr;
	GiNaC::numeric argg;
	GiNaC::ex resultg, resultgRe, resultgIm;
	ostringstream oss,oss2;
/*
	Form measures precision in bits, GiNaC in decimal digits
*/
	int numdigits = (AC.DefaultPrecision-AC.MaxWeight)*log10(2.0);
	GiNaC::Digits = numdigits;
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
	resultgRe = GiNaC::real_part(resultg);
	resultgIm = GiNaC::imag_part(resultg);
/*
	Convert back to mpf_t
*/
	oss << resultgRe;
	mpf_set_str(resultRe, oss.str().c_str(), 10);
	oss2 << resultgIm;
	mpf_set_str(resultIm, oss2.str().c_str(), 10);
	return(0);
}
/*
  	#] CalculateLin :
  	#[ CalculateHpl :
*/
int CalculateHpl(mpf_t resultRe, mpf_t resultIm, int *indexes, int depth, mpf_t arg) {
	int i;
	GiNaC::lst indexesg;
	char *str = nullptr;
	GiNaC::numeric argg;
	GiNaC::ex resultg,resultgRe, resultgIm;
	ostringstream oss,oss2;
/*
	Form measures precision in bits, GiNaC in decimal digits
*/
	int numdigits = (AC.DefaultPrecision-AC.MaxWeight)*log10(2.0);
	GiNaC::Digits = numdigits;
/*
	Trailing zeroes lead to factors of hpl_(0;x) = ln_(x)
	So we only consider x != 0 in this case
*/
	if ( (indexes[depth-1] == 0) && (mpf_cmp_ui(arg,0L) == 0) ) {
		return(-1);
	}
/*
    Leading ones lead to factors of hpl_(1;x) = ln_(1-x) 
	So we only consider x =! 1 in this case
*/
	if ( (indexes[0] == 1) && (mpf_cmp_ui(arg,1L) == 0) ) {
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
	resultgRe = GiNaC::real_part(resultg);
	resultgIm = GiNaC::imag_part(resultg);
/*
	Convert back to mpf_t
*/
	oss << resultgRe;
	mpf_set_str(resultRe, oss.str().c_str(), 10);
	oss2 << resultgIm;
	mpf_set_str(resultIm, oss2.str().c_str(), 10);
	return(0);
}
/*
  	#] CalculateHpl :
  	#[ CalculateMpl :
*/
int CalculateMpl(mpf_t resultRe, mpf_t resultIm, int *indexes, int depth, mpf_t *arg) {
	int i;
	GiNaC::lst indexesg,argListg;
	char *str = nullptr;
	GiNaC::numeric argg;
	GiNaC::ex resultg,resultgRe, resultgIm;
	ostringstream oss,oss2;
/*
	Form measures precision in bits, GiNaC in decimal digits
*/
	int numdigits = (AC.DefaultPrecision-AC.MaxWeight)*log10(2.0);
	GiNaC::Digits = numdigits;
/*
	Place the indexes and arguments in a ginac list
*/
	for( i=0; i < depth; i++ ) {
		indexesg.append(indexes[i]);
		gmp_asprintf(&str, "%.Ff", arg[i]);
    	argg = GiNaC::numeric(str);
		free(str);
		argListg.append(argg);
	}
/*
	Convert to ginac numeric and evaluate
*/
	resultg = GiNaC::evalf(GiNaC::Li(indexesg, argListg));
	resultgRe = GiNaC::real_part(resultg);
	resultgIm = GiNaC::imag_part(resultg);
/*
	Convert back to mpf_t
*/
	oss << resultgRe;
	mpf_set_str(resultRe, oss.str().c_str(), 10);
	oss2 << resultgIm;
	mpf_set_str(resultIm, oss2.str().c_str(), 10);
	return(0);
}
/*
  	#] CalculateMpl :
  	#[ EvaluatePolylog : 
*/

int EvaluatePolylog(PHEAD WORD *term, WORD level, WORD par) {
	WORD *t, *tstop, *tt, *newterm, i, depth, first = 1;
	WORD *oldworkpointer = AT.WorkPointer, nsize, *indexes;
	int retval;
	mpf_t tmpRe,tmpIm, tmp;
	mpf_init(tmpRe); mpf_init(tmpIm); mpf_init(tmp);

	tstop = term + *term; tstop -= ABS(tstop[-1]);
	if ( AT.WorkPointer < term+*term ) AT.WorkPointer = term + *term;

/*
	Step 1: locate a LINFUNCTION, HPLFUNCTION or MPLFUNCTION
*/
	t = term+1;
	while ( t < tstop ) {
		if ( (*t == par) || ( ( par == ALLPOLYLOGFUNCTIONS ) && 
				( *t == LINFUNCTION || *t == HPLFUNCTION || *t == MPLFUNCTION ) ) ) {
			indexes = AT.WorkPointer;
			if (GetLinArgument(indexes,&depth,mpftab1,t) != 0) {
				MesPrint("Error: LINFUNCTION with illegal argument(s)");
				goto nextfun;
			}
/*
	Step 2: evaluate
	The real part is accumulated in aux4, the imaginary part in aux5.
*/
			if ( first ) {
				if ( *t == LINFUNCTION ) {
					if (CalculateLin(aux4,aux5,indexes[0],*mpftab1) != 0) goto nextfun;
				}
				else if ( *t == HPLFUNCTION ) {
					if (CalculateHpl(aux4,aux5,indexes,depth,*mpftab1) != 0) goto nextfun;
				}
				else if ( *t == MPLFUNCTION ) {
					if (CalculateMpl(aux4,aux5,indexes,depth,mpftab1) != 0) goto nextfun;
				}
				first = 0;
			}
			else {
				if ( *t == LINFUNCTION ) {
					if (CalculateLin(auxjm,auxjjm,indexes[0],*mpftab1) != 0) goto nextfun;
				}
				else if ( *t == HPLFUNCTION ) {
					if (CalculateHpl(auxjm,auxjjm,indexes,depth,*mpftab1) != 0) goto nextfun;
				}
				else if ( *t == MPLFUNCTION ) {
					if (CalculateMpl(auxjm,auxjjm,indexes,depth,mpftab1) != 0) goto nextfun;
				}
/*
	We multiply two complex numbers:
	(aux4+aux5*i_) * (auxjm+auxjjm*i_) 
			= (aux4*auxjm-aux5*auxjjm) + (aux4*auxjjm+aux5*auxjm)*i_
*/
				// Real part
				mpf_mul(tmpRe,aux4,auxjm);
				mpf_mul(tmp,aux5,auxjjm);
				mpf_sub(tmpRe,tmpRe,tmp);
				// Imaginary part
				mpf_mul(tmpIm,aux4,auxjjm);
				mpf_mul(tmp,aux5,auxjm);
				mpf_add(tmpIm,tmpIm,tmp);
				// We set the results again in aux4 and aux5
				mpf_set(aux4,tmpRe);
				mpf_set(aux5,tmpIm);
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
		mpf_neg(aux5,aux5);
		nsize = -nsize;
	}
	if ( nsize == 3 ) {
		if ( tstop[0] != 1 ) {
			mpf_mul_ui(aux4,aux4,(ULONG)((UWORD)tstop[0]));
			mpf_mul_ui(aux5,aux5,(ULONG)((UWORD)tstop[0]));
		}
		if ( tstop[1] != 1 ) {
			mpf_div_ui(aux4,aux4,(ULONG)((UWORD)tstop[1]));
			mpf_div_ui(aux5,aux5,(ULONG)((UWORD)tstop[1]));
		}
	}
	else {
		RatToFloat(auxjjm,(UWORD *)tstop,nsize);
		mpf_mul(aux4,aux4,auxjjm);
		mpf_mul(aux5,aux5,auxjjm);
	}
/*
	Now we have to locate possible other float_ functions.
*/
	t = term+1;
	while ( t < tstop ) {
		if ( *t == FLOATFUN ) {
			UnpackFloat(auxjjm,t);
			mpf_mul(aux4,aux4,auxjjm);
			mpf_mul(aux5,aux5,auxjjm);
		}
		t += t[1];
	}

/*
	Now we should compose the new term in the WorkSpace.
*/
	// Real part
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

	// Imaginary part
	if ( mpf_cmp_ui(aux5,0L) != 0 ) {
	t = term+1;
	newterm = AT.WorkPointer;
	tt = newterm+1;
	while ( t < tstop ) {
		if ( *t == 0 || *t == FLOATFUN ) t += t[1];
		else {
			i = t[1]; NCOPY(tt,t,i);
		}
	}
	*tt++ = SYMBOL; *tt++ = 4; *tt++ = ISYMBOL; *tt++ = 1;
	PackFloat(tt,aux5);
	tt += tt[1];
	*tt++ = 1; *tt++ = 1; *tt++ = 3;
	*newterm = tt-newterm;
	AT.WorkPointer = tt;
	retval = Generator(BHEAD newterm,level);
	}
	AT.WorkPointer = oldworkpointer;
	mpf_clear(tmpRe); mpf_clear(tmpIm); mpf_clear(tmp);
	return(retval);
}
/*
  	#] EvaluatePolylog :
*/
