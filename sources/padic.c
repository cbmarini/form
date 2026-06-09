/** @file padic.c
 *
 *  This file contains the p-adic runtime integration.
 *  The implementation follows the same integration points as float.c:
 *  - a dedicated internal function `padic_` used as coefficient carrier,
 *  - statement support (`ToPadic`) in compiler/executor,
 *  - normalization/sorting helper routines for coefficient arithmetic,
 *  - print support with `Format padicprecision` and `Format padicformat`.
 *
 *  The numerical backend is FLINT's padic module.
 */
/* #[ License : */
/*
 *   Copyright (C) 1984-2026 J.A.M. Vermaseren
 *   When using this file you are requested to refer to the publication
 *   J.A.M.Vermaseren "New features of FORM" math-ph/0010025
 *   This is considered a matter of courtesy as the development was paid
 *   for by FOM the Dutch physics granting agency and we would like to
 *   be able to track its scientific use to convince FOM of its value
 *   for the community.
 *
 *   This file is part of FORM.
 *
 *   FORM is free software: you can redistribute it and/or modify it under the
 *   terms of the GNU General Public License as published by the Free Software
 *   Foundation, either version 3 of the License, or (at your option) any later
 *   version.
 *
 *   FORM is distributed in the hope that it will be useful, but WITHOUT ANY
 *   WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 *   FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 *   details.
 *
 *   You should have received a copy of the GNU General Public License along
 *   with FORM.  If not, see <http://www.gnu.org/licenses/>.
 */
/* #] License : */
/*
  	#[ Includes: padic.c
*/
#include "form3.h"
#include <stdio.h>
#include <string.h>

/*
 * FLINT's padic mpq conversion API requires gmp.h to be included before
 * flint headers.
 */
#include <gmp.h>
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpq.h>
#include <flint/padic.h>

/*
	FORM keeps a single active p-adic context per run:
	- prime p
	- precision N
	The context is configured by #StartPadic / #EndPadic.

	Note: the internal `padic_` coefficient carrier does not store p itself.
	The active `PadicPrime`/`PadicContext` is assumed when unpacking and
	operating on p-adic coefficients.
*/
#define PadicActive            AC.activePadic
#define PadicPrime             AC.activePadicPrime
#define PadicPrecision         AC.activePadicPrecision
#define ActivePadicContext     AC.activePadicContext
#define PadicContext           ((padic_ctx_struct *)(ActivePadicContext))

/*
	Thread-local aux storage. The pointer is stored in AT.padic_aux_.
*/
typedef struct PADIC_AUX_ {
	/*
		All p-adic operations are done through this per-thread scratch space.
		This avoids repeated init/clear churn while sorting/normalizing and keeps
		GMP/FLINT temporaries thread-local.
	*/
	padic_t p1;
	padic_t p2;
	padic_t p3;
	/* Used by PackPadic() to build a canonical (reduced) representation. */
	padic_t p4;
	/* Scratch rational used by FORM <-> padic/fmpq conversions. */
	fmpq_t q1;
	/* Scratch integer used by list-format printing. */
	mpz_t z2;
} PADIC_AUX;

/*
	AT.padic_aux_ is allocated by StartPadicSystem() for all threads.
	When p-adics are not active this pointer is 0.
*/
#define PadicAux ((PADIC_AUX *)(AT.padic_aux_))
#define paux1 (PadicAux->p1)
#define paux2 (PadicAux->p2)
#define paux3 (PadicAux->p3)
#define paux4 (PadicAux->p4)
#define pauxq1 (PadicAux->q1)
#define pauxz2 (PadicAux->z2)

/*
 	#] Includes : 
  	#[ Helpers :
 		#[ InitPadicAux :
*/
static void InitPadicAux(PADIC_AUX *aux, LONG prec)
{
	padic_init2(aux->p1, (slong)prec);
	padic_init2(aux->p2, (slong)prec);
	padic_init2(aux->p3, (slong)prec);
	padic_init2(aux->p4, (slong)prec);
	fmpq_init(aux->q1);
	mpz_init(aux->z2);
}
/*
 		#] InitPadicAux :
 		#[ ClearSinglePadicAux :
*/
static void ClearSinglePadicAux(PADIC_AUX *aux)
{
	padic_clear(aux->p1);
	padic_clear(aux->p2);
	padic_clear(aux->p3);
	padic_clear(aux->p4);
	fmpq_clear(aux->q1);
	mpz_clear(aux->z2);
}
/*
 		#] ClearSinglePadicAux :
 		#[ AllocatePadicAux :
*/
static void AllocatePadicAux(void)
{
#ifdef WITHPTHREADS
	int id, totnum;
	totnum = AM.totalnumberofthreads;
#ifdef WITHSORTBOTS
	totnum = MaX(2 * AM.totalnumberofthreads - 3, AM.totalnumberofthreads);
#endif
	for ( id = 0; id < totnum; id++ ) {
		PADIC_AUX *aux;
		if ( AB[id]->T.padic_aux_ ) {
			ClearSinglePadicAux((PADIC_AUX *)(AB[id]->T.padic_aux_));
			M_free(AB[id]->T.padic_aux_,"AB[id]->T.padic_aux_");
			AB[id]->T.padic_aux_ = 0;
		}
		aux = Malloc1(sizeof(PADIC_AUX),"AB[id]->T.padic_aux_");
		InitPadicAux(aux,PadicPrecision);
		AB[id]->T.padic_aux_ = (void *)aux;
	}
#else
	PADIC_AUX *aux;
	/*
		Single-thread build: AT.padic_aux_ is the only instance.
	*/
	if ( AT.padic_aux_ ) {
		ClearSinglePadicAux((PADIC_AUX *)(AT.padic_aux_));
		M_free(AT.padic_aux_,"AT.padic_aux_");
		AT.padic_aux_ = 0;
	}
	aux = Malloc1(sizeof(PADIC_AUX),"AT.padic_aux_");
	InitPadicAux(aux,PadicPrecision);
	AT.padic_aux_ = (void *)aux;
#endif
}
/*
 		#] AllocatePadicAux :
 		#[ ClearPadicAux :
*/
static void ClearPadicAux(void)
{
#ifdef WITHPTHREADS
	int id, totnum;
	totnum = AM.totalnumberofthreads;
#ifdef WITHSORTBOTS
	totnum = MaX(2 * AM.totalnumberofthreads - 3, AM.totalnumberofthreads);
#endif
	for ( id = 0; id < totnum; id++ ) {
		if ( AB[id]->T.padic_aux_ ) {
			ClearSinglePadicAux((PADIC_AUX *)(AB[id]->T.padic_aux_));
			M_free(AB[id]->T.padic_aux_,"AB[id]->T.padic_aux_");
			AB[id]->T.padic_aux_ = 0;
		}
	}
#else
	if ( AT.padic_aux_ ) {
		ClearSinglePadicAux((PADIC_AUX *)(AT.padic_aux_));
		M_free(AT.padic_aux_,"AT.padic_aux_");
		AT.padic_aux_ = 0;
	}
#endif
}
/*
 		#] ClearPadicAux :
 		#[ PadicIsPrime :
*/
int PadicIsPrime(LONG p)
{
	int prime;
	fmpz_t z;
	if ( p <= 1 ) return(0);
	fmpz_init_set_si(z,(slong)p);
	prime = fmpz_is_prime(z);
	fmpz_clear(z);
	return(prime == 1);
}
/*
 		#] PadicIsPrime :
 		#[ StartPadicSystem :

	Initializes (or reinitializes) the single global p-adic context.
	Called by #StartPadic from the preprocessor.

	This function:
	- clears the previous context if active,
	- initializes FLINT's padic_ctx_struct for (p,N),
	- allocates per-thread scratch objects (AT.padic_aux_ / AB[id]->T.padic_aux_).
*/
int StartPadicSystem(LONG p, LONG N)
{
	fmpz_t prime;
	if ( p <= 1 || N <= 0 ) return(1);
	if ( PadicActive ) { // Clear the previous padic system
		ClearPadicSystem();
	}
	PadicPrime = p;
	PadicPrecision = N;
	ActivePadicContext = Malloc1(sizeof(padic_ctx_struct),"PadicContext");
	fmpz_init_set_si(prime,(slong)p);
	padic_ctx_init(PadicContext,prime,0,(slong)N,PADIC_SERIES);
	fmpz_clear(prime);
	AllocatePadicAux();
	PadicActive = 1;
	return(0);
}
/*
 		#] StartPadicSystem :
 		#[ ClearPadicSystem :

	Releases all runtime state associated with p-adic arithmetic, including:
	- per-thread aux buffers,
	- the FLINT context,
	- cached print buffer AO.padicspace.
*/
void ClearPadicSystem(void)
{
	ClearPadicAux();
	if ( PadicActive ) {
		padic_ctx_clear(PadicContext);
	}
	if ( ActivePadicContext ) {
		M_free(ActivePadicContext,"PadicContext");
		ActivePadicContext = 0;
	}
	PadicActive = 0;
	PadicPrime = 0;
	PadicPrecision = 0;
	if ( AO.padicspace ) {
		M_free(AO.padicspace,"padicspace");
		AO.padicspace = 0;
		AO.padicsize = 0;
	}
}
/*
 		#] ClearPadicSystem :
  	#] Helpers :
  	#[ Internal p-adic function format :
 		#[ Explanations :

	p-adic coefficients are stored in a special function padic_ whose 
	arguments contain the relevant p-adic information from the FLINT library:
	    padic_(v, N, u)
	where
	- v is the p-adic valuation,
	- N is the precision stored with this coefficient,
	- u is the FLINT unit.

	The prime p is not stored in the function. It is supplied by the active 
	p-adic context configured through #StartPadic / #EndPadic. The unit
	u is stored in the numerator of a normal FORM integer with denominator 1. 
	Note: we could store p in the denominator as p cannot divide u.

		#] Explanations :
 		#[ TestPadic :

	Checks whether fun has the shape of a legal padic_ with its arguments:
	    padic_(v,N,u)

	Return value:
	- 1 if fun is a syntactically valid padic_ function,
	- 0 otherwise.
*/
int TestPadic(WORD *fun)
{
	WORD *f, *fstop;
	WORD nargs, nnum, i;
	LONG N;
	ULONG x;

	f = fun + FUNHEAD;
	fstop = fun + fun[1];
	nargs = 0;
	while ( f < fstop ) {
		nargs++;
		NEXTARG(f);
	}
	if ( nargs != 3 || f != fstop ) return(0);

	f = fun + FUNHEAD;
	/* v: signed LONG argument */
	if ( *f == -SNUMBER ) {
		f += 2;
	}
	else {
		if ( *f != ARGHEAD+6 ) return(0);
		if ( ABS(f[ARGHEAD+5]) != 5 ) return(0);
		if ( f[ARGHEAD+3] != 1 ) return(0);
		if ( f[ARGHEAD+4] != 0 ) return(0);
		f += *f;
	}
	/* N: signed LONG argument; additionally require positive precision. */
	if ( *f == -SNUMBER ) {
		N = (LONG)f[1];
		f += 2;
	}
	else {
		if ( *f != ARGHEAD+6 ) return(0);
		if ( ABS(f[ARGHEAD+5]) != 5 ) return(0);
		if ( f[ARGHEAD+3] != 1 ) return(0);
		if ( f[ARGHEAD+4] != 0 ) return(0);
		x = ((ULONG)(UWORD)f[ARGHEAD+2] << BITSINWORD) + (UWORD)f[ARGHEAD+1];
		N = (f[ARGHEAD+5] < 0) ? -(LONG)x : (LONG)x;
		f += *f;
	}
	if ( N <= 0 ) return(0);

	/* u: FORM integer (-SNUMBER or long integer n/1). */
	if ( *f == -SNUMBER ) {
		f += 2;
	}
	else {
		nnum = (WORD)((ABS(f[*f-1])-1)/2);
		if ( nnum <= 0 ) return(0);
		if ( *f != ARGHEAD + 2*nnum + 2 ) return(0);
		if ( f[ARGHEAD] != 2*nnum + 2 ) return(0);
		/* Denominator must be exactly 1. */
		f += ARGHEAD + nnum + 1;
		if ( f[0] != 1 ) return(0);
		for ( i = 1; i < nnum; i++ ) {
			if ( f[i] != 0 ) return(0);
		}
	}
	return(1);
}
/*
 		#] TestPadic :
 		#[ UnpackPadic :

	Converts the arguments of a padic_ function back into FLINT's padic_t.
	We shouldn't come here if the p-adic system is turned off. 

	Return value:
	- 0  on success,
	- -1 if fun is not a valid internal p-adic function.
*/
static int UnpackPadic(padic_t out, WORD *fun)
{
	WORD *f;
	ULONG x;

	if ( !PadicActive ) {
		MLOCK(ErrorMessageLock);
		MesPrint("Illegal attempt at using a padic_ function without proper startup.");
		MesPrint("Please use %#StartPadic <p>,<N> first.");
		MUNLOCK(ErrorMessageLock);
		Terminate(-1);
	}

	if ( TestPadic(fun) == 0 ) return(-1);

	/*
		Read the argument triplet in order: v, N, u.
	*/
	f = fun + FUNHEAD;

	/*
		TestPadic() already validated the arguments of padic_, 
		so we can safely decode v, N and u here.
	*/
	if ( *f == -SNUMBER ) {
		padic_val(out) = (slong)f[1];
		f += 2;
	}
	else {
		x = ((ULONG)(UWORD)f[ARGHEAD+2] << BITSINWORD) + (UWORD)f[ARGHEAD+1];
		padic_val(out) = (f[ARGHEAD+5] < 0) ? -(slong)x : (slong)x;
		f += *f;
	}

	if ( *f == -SNUMBER ) {
		padic_prec(out) = (slong)f[1];
		f += 2;
	}
	else {
		x = ((ULONG)(UWORD)f[ARGHEAD+2] << BITSINWORD) + (UWORD)f[ARGHEAD+1];
		padic_prec(out) = (f[ARGHEAD+5] < 0) ? -(slong)x : (slong)x;
		f += *f;
	}

	/*
		Decode the unit u.
	*/
	if ( *f == -SNUMBER ) {
		fmpz_set_si(padic_unit(out),(slong)(f[1]));
	}
	else {
		WORD size, nnum;
		size = f[*f-1];
		nnum = (WORD)((ABS(size)-1)/2);
		flint_fmpz_set_form(padic_unit(out),(UWORD *)(f+ARGHEAD+1),
			(size < 0) ? -nnum : nnum);
	}

	padic_reduce(out,PadicContext);
	return(0);
}
/*
 		#] UnpackPadic :
 		#[ PackPadic :

	Converts a FLINT's padic_t into FORM's internal padic notation:
	    padic_(v,N,u)

	Return value:
	- the number of WORDs written to fun.
*/
static int PackPadic(WORD *fun, padic_t in)
{
	WORD *t;
	LONG v, N, small;
	ULONG x;
	WORD *arg;
	WORD nnum, nabs, i;
	GETIDENTITY

	/*
		Normalize first so the (v,N,u) triplet is canonical.
	*/
	padic_set(paux4,in,PadicContext);
	padic_reduce(paux4,PadicContext);
	v = (LONG)padic_val(paux4);
	N = (LONG)padic_prec(paux4);

	/*
		We now fill the function with the three arguments:
		valuation v, precision N, and unit u.
	*/
	t = fun;
	*t++ = PADICFUN;
	t++; // fun[1] (function length) is filled at the end.
	FILLFUN(t);

	/*
		Pack valuation v:
		compact -SNUMBER when it fits in WORD, otherwise as a two-word
		long integer argument.
	*/
	if ( v >= WORD_MIN_VALUE && v <= WORD_MAX_VALUE ) {
		*t++ = -SNUMBER;
		*t++ = (WORD)v;
	}
	else {
		x = ( v < 0 ) ? (ULONG)(-v) : (ULONG)v;
		*t++ = ARGHEAD + 6;
		*t++ = 0;
		FILLARG(t);
		*t++ = 6;
		*t++ = (UWORD)x;
		*t++ = (UWORD)(x >> BITSINWORD);
		*t++ = 1;
		*t++ = 0;
		*t++ = (v < 0) ? -5 : 5;
	}
	/*
		Pack precision N.
	*/
	if ( N >= WORD_MIN_VALUE && N <= WORD_MAX_VALUE ) {
		*t++ = -SNUMBER;
		*t++ = (WORD)N;
	}
	else {
		x = ( N < 0 ) ? (ULONG)(-N) : (ULONG)N;
		*t++ = ARGHEAD + 6;
		*t++ = 0;
		FILLARG(t);
		*t++ = 6;
		*t++ = (UWORD)x;
		*t++ = (UWORD)(x >> BITSINWORD);
		*t++ = 1;
		*t++ = 0;
		*t++ = (N < 0) ? -5 : 5;
	}
	/*
		Pack unit u in the numerator of a FORM integer.
		use -SNUMBER for small values, otherwise FORM long-integer.
	*/
	if ( fmpz_fits_si(padic_unit(paux4))
	  && (small = (LONG)fmpz_get_si(padic_unit(paux4)),
	      small >= WORD_MIN_VALUE && small <= WORD_MAX_VALUE) ) {
		*t++ = -SNUMBER;
		*t++ = (WORD)small;
	}
	else {
		arg = t + ARGHEAD + 1;
		// Numerator
		nnum = flint_fmpz_get_form(padic_unit(paux4),arg);
		nabs = ABS(nnum);
		// Denominator 
		arg += nabs;
		*arg++ = 1;
		for ( i = 1; i < nabs; i++ ) *arg++ = 0;
		*arg++ = ( nnum < 0 ) ? -(2*nabs+1) : (2*nabs+1);
		// Arghead etc.
		*t++ = ARGHEAD + 2*nabs + 2;
		*t++ = 0; 
		FILLARG(t);
		*t++ = 2*nabs + 2;
		t = arg;
	}
	fun[1] = t - fun;
	return(fun[1]);
}
/*
 		#] PackPadic :
  	#] Internal p-adic function format :
  	#[ Rekenen :
		#[ FormRatToFmpq :

	Converts the internal FORM rational coefficient encoding (formrat/ratsize)
	to a FLINT rational fmpq.

	Maybe this should be moved to flintinterface.cc?
*/
static void FormRatToFmpq(fmpq_t result, UWORD *formrat, WORD ratsize)
{
	WORD nnum, nden;
	UWORD *num, *den;
	int sign = 1;

	if ( ratsize < 0 ) {
		sign = -1;
		ratsize = -ratsize;
	}

	nnum = nden = (ratsize-1)/2;
	num = formrat;
	den = formrat + nnum;
	/* Remove padding */
	while ( nnum > 0 && num[nnum-1] == 0 ) nnum--;
	while ( nden > 0 && den[nden-1] == 0 ) nden--;

	if ( nnum > 0 ) {
		flint_fmpz_set_form(fmpq_numref(result),num,sign*nnum);
	}
	else {
		fmpz_zero(fmpq_numref(result));
	}
	flint_fmpz_set_form(fmpq_denref(result),den,nden);

	fmpq_canonicalise(result);
}
/*
		#] FormRatToFmpq :
 		#[ MulRatToPadic :

	Multiplies an existing `padic_` coefficient by a FORM rational coefficient.
	This is used by Normalize() when a term contains both a numeric coefficient
	and a `padic_` function.

	Return value:
	- 0  on success, with a non-zero product packed into outfun,
	- 1  if the product is p-adic zero,
	- -1 on runtime/validation failure.
*/
int MulRatToPadic(PHEAD WORD *outfun, WORD *infun, UWORD *formrat, WORD nrat)
{
	if ( !PadicActive ) return(-1);
	if ( UnpackPadic(paux1,infun) ) return(-1);
	FormRatToFmpq(pauxq1,formrat,nrat);
	padic_set_fmpq(paux2,pauxq1,PadicContext);
	padic_mul(paux3,paux1,paux2,PadicContext);
	PackPadic(outfun,paux3);
	if ( padic_is_zero(paux3) ) return(1);
	return(0);
}
/*
 		#] MulRatToPadic :
 		#[ MulPadics :

	Multiplies two padic_ functions (fun1 * fun2) and stores the
	product as a new padic_ function in fun3.

	Return value:
	- 0  on success, with a non-zero product packed into fun3,
	- 1  if the product is p-adic zero,
	- -1 on runtime/validation failure.
*/
int MulPadics(PHEAD WORD *fun3, WORD *fun1, WORD *fun2)
{
	if ( !PadicActive ) return(-1);
	if ( UnpackPadic(paux1,fun1) ) return(-1);
	if ( UnpackPadic(paux2,fun2) ) return(-1);
	padic_mul(paux3,paux1,paux2,PadicContext);
	PackPadic(fun3,paux3);
	if ( padic_is_zero(paux3) ) return(1);
	return(0);
}
/*
 		#] MulPadics :
 		#[ DivPadics :

	Divides two padic_ functions (fun1 / fun2) and stores the
	quotient as a new padic_ function in fun3.
	Division by zero is a fatal runtime error.

	Return value:
	- 0  on success, with a non-zero quotient packed into fun3,
	- 1  if the quotient is p-adic zero,
	- -1 on runtime/validation failure.
*/
int DivPadics(PHEAD WORD *fun3, WORD *fun1, WORD *fun2)
{
	if ( !PadicActive ) return(-1);
	if ( UnpackPadic(paux1,fun1) ) return(-1);
	if ( UnpackPadic(paux2,fun2) ) return(-1);
	if ( padic_is_zero(paux2) ) {
		MLOCK(ErrorMessageLock);
		MesPrint("Division by zero in p-adic arithmetic.");
		MUNLOCK(ErrorMessageLock);
		Terminate(-1);
		return(-1);
	}
	padic_div(paux3,paux1,paux2,PadicContext);
	PackPadic(fun3,paux3);
	if ( padic_is_zero(paux3) ) return(1);
	return(0);
}
/*
 		#] DivPadics :
		#[ InvPadic :

	Computes the inverse 1/fun and stores the result as a new padic_ function
	in outfun. Division by zero is a fatal runtime error.

	Return value:
	- 0  on success, with a non-zero inverse packed into outfun,
	- 1  if the inverse is p-adic zero,
	- -1 on runtime/validation failure.
*/
int InvPadic(PHEAD WORD *outfun, WORD *fun)
{
	if ( !PadicActive ) return(-1);
	if ( UnpackPadic(paux1,fun) ) return(-1);
	if ( padic_is_zero(paux1) ) {
		MLOCK(ErrorMessageLock);
		MesPrint("Division by zero in p-adic arithmetic.");
		MUNLOCK(ErrorMessageLock);
		Terminate(-1);
		return(-1);
	}
	padic_inv(paux3,paux1,PadicContext);
	PackPadic(outfun,paux3);
	if ( padic_is_zero(paux3) ) return(1);
	return(0);
}
/*
		#] InvPadic :
		#[ FmpqToFormRat :

	Converts a FLINT rational to FORM's internal rational coefficient encoding.

	If out is 0, this routine only computes the signed size code in *nratout.
	Returns -1 when the result does not fit in FORM's encoding bounds.
*/
static int FmpqToFormRat(UWORD *out, WORD *nratout, fmpq_t q)
{
	int sign;
	size_t nnum, nden, n, i;
	fmpz_t z;

	sign = fmpq_sgn(q);
	if ( sign == 0 ) {
		*nratout = 3;
		if ( out != 0 ) {
			out[0] = 0;
			out[1] = 1;
			out[2] = 3;
		}
		return(0);
	}

	nnum = fmpz_is_zero(fmpq_numref(q)) ? 0 :
		(size_t)((fmpz_sizeinbase(fmpq_numref(q),2) + BITSINWORD - 1) / BITSINWORD);
	nden = fmpz_is_zero(fmpq_denref(q)) ? 0 :
		(size_t)((fmpz_sizeinbase(fmpq_denref(q),2) + BITSINWORD - 1) / BITSINWORD);
	if ( nden == 0 ) return(-1);
	n = ( nnum > nden ) ? nnum : nden;
	if ( n > (size_t)((WORD_MAX_VALUE-1)/2) ) return(-1);

	*nratout = (WORD)(2*n + 1);
	if ( sign < 0 ) *nratout = -*nratout;
	if ( out == 0 ) return(0);

	for ( i = 0; i < 2*n; i++ ) out[i] = 0;

	fmpz_init(z);
	if ( nnum > 0 ) {
		fmpz_abs(z,fmpq_numref(q));
		flint_fmpz_get_form(z,(WORD *)out);
	}
	if ( nden > 0 ) {
		fmpz_set(z,fmpq_denref(q));
		flint_fmpz_get_form(z,(WORD *)(out+n));
	}
	fmpz_clear(z);
	out[2*n] = (UWORD)ABS(*nratout);
	return(0);
}
/*
		#] FmpqToFormRat :
		#[ PadicReconstructToFmpq :

		Reconstructs a small rational from a reduced p-adic value using FLINT's
		rational reconstruction and applies the p-adic valuation afterwards.

		This is a "best effort" conversion used by PadicToRat(): if reconstruction
		fails (not unique for the current modulus), the caller can fall back to a
		canonical lift via padic_get_fmpq().
 */
static int PadicReconstructToFmpq(fmpq_t out, padic_t in)
{
	fmpz_t residue, modulus, ppower;
	fmpq_t recon;
	slong v, N, modexp;
	int ok = 0;

	fmpz_init(residue);
	fmpz_init(modulus);
	fmpz_init(ppower);
	fmpq_init(recon);

	v = padic_get_val(in);
	N = padic_get_prec(in);
	modexp = N - v;
	if ( modexp <= 0 ) goto ClearAndReturn;

	fmpz_pow_ui(modulus,PadicContext->p,(ulong)modexp);
	fmpz_mod(residue,padic_unit(in),modulus);
	ok = fmpq_reconstruct_fmpz(recon,residue,modulus);
	if ( !ok ) goto ClearAndReturn;

	if ( v > 0 ) {
		fmpz_pow_ui(ppower,PadicContext->p,(ulong)v);
		fmpq_mul_fmpz(recon,recon,ppower);
	}
	else if ( v < 0 ) {
		fmpz_pow_ui(ppower,PadicContext->p,(ulong)(-v));
		fmpq_div_fmpz(recon,recon,ppower);
	}

	fmpq_set(out,recon);

ClearAndReturn:
	fmpq_clear(recon);
	fmpz_clear(ppower);
	fmpz_clear(modulus);
	fmpz_clear(residue);
	return(ok ? 0 : -1);
}
/*
		#] PadicReconstructToFmpq :
  	#] Rekenen :
  	#[ Printing :
 		#[ EnsurePadicPrintBuffer :
*/
static int EnsurePadicPrintBuffer(size_t need)
{
	if ( AO.padicspace == 0 || AO.padicsize <= (LONG)need ) {
		if ( AO.padicspace ) M_free(AO.padicspace,"padicspace");
		AO.padicsize = (LONG)need + 32;
		AO.padicspace = (UBYTE *)Malloc1((size_t)AO.padicsize,"padicspace");
	}
	return(AO.padicspace == 0 ? -1 : 0);
}
/*
 		#] EnsurePadicPrintBuffer :
 		#[ CountULongDigits :
*/
static size_t CountULongDigits(unsigned long x)
{
	size_t n = 1;
	while ( x >= 10 ) {
		x /= 10;
		n++;
	}
	return(n);
}
/*
 		#] CountULongDigits :
 		#[ PrintPadicList :
*/
static int PrintPadicList(PADIC_AUX *aux, padic_t in)
{
	char prefix[96];
	char *out;
	slong v, prec, i, coeff_count;
	size_t n, prefix_len, total_len;
	mpz_t work, rem, pz;

	v = padic_val(in);
	prec = padic_prec(in);
	coeff_count = prec - v;

	/*
		The textual list is a plain coefficient list without Big-O tail.
		When the truncated value is effectively zero, print a canonical zero.
	*/
	if ( fmpz_is_zero(padic_unit(in)) || coeff_count <= 0 ) {
		n = (size_t)snprintf(prefix,sizeof(prefix),"padic[%ld,0,%ld,{0}]",
			(long)PadicPrime,(long)PadicPrecision);
		if ( n >= sizeof(prefix) ) return(0);
		if ( EnsurePadicPrintBuffer(n) ) return(0);
		memcpy((char *)AO.padicspace,prefix,n+1);
		return((int)n);
	}

	fmpz_get_mpz(aux->z2,padic_unit(in));
	mpz_init_set(work,aux->z2);
	mpz_init(rem);
	mpz_init_set_ui(pz,(unsigned long)PadicPrime);

	prefix_len = (size_t)snprintf(prefix,sizeof(prefix),"padic[%ld,%ld,%ld,{",
		(long)PadicPrime,(long)v,(long)PadicPrecision);
	if ( prefix_len >= sizeof(prefix) ) {
		mpz_clear(pz);
		mpz_clear(rem);
		mpz_clear(work);
		return(0);
	}

	total_len = prefix_len + 2; /* "}]" */
	for ( i = 0; i < coeff_count; i++ ) {
		unsigned long coeff;
		mpz_fdiv_qr(work,rem,work,pz);
		coeff = mpz_get_ui(rem);
		total_len += CountULongDigits(coeff);
		if ( i + 1 < coeff_count ) total_len++;
	}

	if ( EnsurePadicPrintBuffer(total_len) ) {
		mpz_clear(pz);
		mpz_clear(rem);
		mpz_clear(work);
		return(0);
	}

	out = (char *)AO.padicspace;
	memcpy(out,prefix,prefix_len);
	n = prefix_len;

	mpz_set(work,aux->z2);
	for ( i = 0; i < coeff_count; i++ ) {
		unsigned long coeff;
		size_t wrote;
		mpz_fdiv_qr(work,rem,work,pz);
		coeff = mpz_get_ui(rem);
		wrote = (size_t)snprintf(out + n,(size_t)(AO.padicsize - (LONG)n),"%lu",coeff);
		n += wrote;
		if ( i + 1 < coeff_count ) out[n++] = ',';
	}
	out[n++] = '}';
	out[n++] = ']';
	out[n] = 0;

	mpz_clear(pz);
	mpz_clear(rem);
	mpz_clear(work);
	return((int)n);
}
/*
 		#] PrintPadicList :
 		#[ PrintPadic :

	Formats a padic_ coefficient for printing.

	Two output formats are supported:
	- series format (default): FLINT's p-adic series text,
	- list format:            padic[p,v,N,{q1,q2,...}] where:
	  * p is the active prime,
	  * v is the valuation of the printed value,
	  * N is the active context precision,
	  * q_i are the series coefficients in ascending powers of p.

	The resulting C string is stored in AO.padicspace and the return value is
	the string length. FORM's print backend reads AO.padicspace after this call.
	Series mode keeps surrounding parentheses; list mode prints plain padic[...].
*/
int PrintPadic(WORD *fun,int numdigits)
{
	GETIDENTITY
	PADIC_AUX *aux;
	char *flint_string;
	size_t n;
	int digits = (int)PadicPrecision;
	int mode = AO.PadicFormat;

	if ( !PadicActive ) return(0);
	aux = PadicAux;
	if ( UnpackPadic(aux->p1,fun) ) return(0);

	if ( numdigits > 0 && numdigits < digits ) digits = numdigits;

	if ( digits == (int)PadicPrecision ) {
		if ( mode == PADICPRINTLIST ) return(PrintPadicList(aux,aux->p1));
		flint_string = padic_get_str(0,aux->p1,PadicContext);
	}
	else {
		padic_ctx_t short_ctx;
		padic_t short_x;
		padic_ctx_init(short_ctx,PadicContext->p,0,(slong)digits,PADIC_SERIES);
		padic_init2(short_x,digits);
		padic_get_fmpq(aux->q1,aux->p1,PadicContext);
		padic_set_fmpq(short_x,aux->q1,short_ctx);
		if ( mode == PADICPRINTLIST ) {
			int outlen = PrintPadicList(aux,short_x);
			padic_clear(short_x);
			padic_ctx_clear(short_ctx);
			return(outlen);
		}
		flint_string = padic_get_str(0,short_x,short_ctx);
		padic_clear(short_x);
		padic_ctx_clear(short_ctx);
	}
	if ( flint_string == 0 ) return(0);

	n = strlen(flint_string);
	/* When truncating digits, append a Big-O tail unless FLINT already did. */
	if ( digits < (int)PadicPrecision && strstr(flint_string,"O(") == 0 ) {
		char tail[64];
		size_t tail_len;
		size_t new_len;
		snprintf(tail,sizeof(tail)," + O(%ld^%d)",(long)PadicPrime,digits);
		tail_len = strlen(tail);
		new_len = n + tail_len + 2;
		if ( EnsurePadicPrintBuffer(new_len) ) { flint_free(flint_string); return(0); }
		((char *)AO.padicspace)[0] = '(';
		memcpy((char *)AO.padicspace + 1, flint_string, n);
		memcpy((char *)AO.padicspace + 1 + n, tail, tail_len);
		((char *)AO.padicspace)[1 + n + tail_len] = ')';
		((char *)AO.padicspace)[new_len] = 0;
		flint_free(flint_string);
		return((int)new_len);
	}

	n += 2;
	if ( EnsurePadicPrintBuffer(n) ) { flint_free(flint_string); return(0); }
	((char *)AO.padicspace)[0] = '(';
	memcpy((char *)AO.padicspace + 1, flint_string, n - 2);
	((char *)AO.padicspace)[n - 1] = ')';
	((char *)AO.padicspace)[n] = 0;
	flint_free(flint_string);
	return((int)n);
}
/*
 		#] PrintPadic :
  	#] Printing :
  	#[ Compiler/runtime statements :
 		#[ CoToPadic :

	Compiler front-end for the `ToPadic;` statement.
	This only validates syntax and records the action; the actual conversion is
	performed on terms at execution time by ToPadic().
*/
int CoToPadic(UBYTE *s)
{
	if ( !PadicActive ) {
		MesPrint("&Illegal attempt to convert to padic_ without activating p-adic numbers.");
		MesPrint("&Forgotten %#startpadic instruction?");
		return(1);
	}
	while ( *s == ' ' || *s == ',' || *s == '\t' ) s++;
	if ( *s ) {
		MesPrint("&Illegal argument(s) in ToPadic statement: '%s'",s);
		return(1);
	}
	Add2Com(TYPETOPADIC);
	return(0);
}
/*
 		#] CoToPadic :
 		#[ CoPadicToRat :

	Compiler front-end for the `PadicToRat;` statement.
	The runtime action is implemented by PadicToRat().
*/
int CoPadicToRat(UBYTE *s)
{
	if ( !PadicActive ) {
		MesPrint("&Illegal attempt to convert from padic_ without activating p-adic numbers.");
		MesPrint("&Forgotten %#startpadic instruction?");
		return(1);
	}
	while ( *s == ' ' || *s == ',' || *s == '\t' ) s++;
	if ( *s ) {
		MesPrint("&Illegal argument(s) in PadicToRat statement: '%s'",s);
		return(1);
	}
	Add2Com(TYPEPADICTORAT);
	return(0);
}
/*
 		#] CoPadicToRat :
 		#[ ToPadic :

	Runtime implementation of `ToPadic;`.

	This replaces the current rational coefficient by an explicit padic_ function
	and sets the rational coefficient to 1/1. The sign is absorbed into the
	p-adic value. If the coefficient is already 1/1 and the term already ends
	with a proper padic_ function we are done immediately.
*/
int ToPadic(PHEAD WORD *term, WORD level)
{
	GETBIDENTITY
	WORD *t, *tstop, nsize, ncoef;

	if ( !PadicActive ) return(1);

	t = term + *term;
	ncoef = t[-1];          
	nsize = ABS(ncoef);     
	tstop = t - nsize;      

	/* If there is already a proper padic_, with rational coefficient 1/1 we are done. */
	if ( ncoef == 3 && t[-2] == 1 && t[-3] == 1 ) {
		t = term + 1;
		while ( t < tstop ) {
			if ( *t == PADICFUN && (t+t[1] == tstop) && TestPadic(t) ) {
				return(Generator(BHEAD term,level));
			}
			t += t[1]; /* advance to the next function in the term */
		}
	}

	FormRatToFmpq(pauxq1,(UWORD *)tstop,ncoef);
	padic_set_fmpq(paux1,pauxq1,PadicContext);
	// Overwrite the old rational coefficient of the term by padic_(v,N,u) 
	// and append rational coefficient 1/1.
	PackPadic(tstop,paux1);
	tstop += tstop[1]; /* advance past the newly written padic_ */
	*tstop++ = 1;
	*tstop++ = 1;
	*tstop++ = 3;
	*term = tstop - term;
	AT.WorkPointer = tstop;
	return(Generator(BHEAD term,level));
}
/*
 		#] ToPadic :
 		#[ PadicToRat :

	Runtime implementation of `PadicToRat;`.

	This finds a terminal `padic_` coefficient record and converts it back to a
	FORM rational coefficient. The conversion first tries rational reconstruction
	(to keep results small) and otherwise falls back to a canonical lift.
*/
int PadicToRat(PHEAD WORD *term, WORD level)
{
	GETBIDENTITY
	WORD *tstop, *t, *stop, nsize, nsign, ncoef;

	if ( !PadicActive ) return(1);

	tstop = term + *term;
	nsize = ABS(tstop[-1]);
	nsign = tstop[-1] < 0 ? -1 : 1;
	tstop -= nsize;
	t = term + 1;
	/*
		The term must end in a single proper padic_ function, followed by a unit
		coefficient 1/1 (the sign is carried separately in tstop[-1]).
	*/
	while ( t < tstop ) {
		if ( *t == PADICFUN && t + t[1] == tstop && TestPadic(t)
		 && nsize == 3 && tstop[0] == 1 && tstop[1] == 1 ) break;
		t += t[1];
	}
	if ( t < tstop ) {
		if ( UnpackPadic(paux1,t) ) return(1);
		if ( padic_is_zero(paux1) ) return(0);
		if ( PadicReconstructToFmpq(pauxq1,paux1) ) {
			/*
				No unique small reconstruction exists for the current precision.
				Fall back to FLINT's canonical lift.
			*/
			padic_get_fmpq(pauxq1,paux1,PadicContext);
		}
		if ( FmpqToFormRat(0,&ncoef,pauxq1) ) goto RatFailure;
		stop = (WORD *)(((UBYTE *)term) + AM.MaxTer);
		if ( t + ABS(ncoef) > stop ) {
			MLOCK(ErrorMessageLock);
			MesPrint("Term too complex after p-adic to rational conversion. MaxTermSize = %10l",
				AM.MaxTer/sizeof(WORD));
			MUNLOCK(ErrorMessageLock);
			Terminate(-1);
			return(1);
		}
		if ( FmpqToFormRat((UWORD *)t,&ncoef,pauxq1) ) goto RatFailure;
		if ( t[0] == 0 && t[1] == 1 && ncoef == 3 ) return(0);
		t += ABS(ncoef);
		t[-1] = ncoef*nsign;
		*term = t - term;
	}
	return(Generator(BHEAD term,level));

RatFailure:
	MLOCK(ErrorMessageLock);
	MesPrint("Failed to convert p-adic coefficient to rational.");
	MUNLOCK(ErrorMessageLock);
	Terminate(-1);
	return(1);
}
/*
 		#] PadicToRat :
  	#] Compiler/runtime statements :
  	#[ Sorting :
 		#[ AddWithPadic :

	Sort helper used when two otherwise identical terms differ only in their
	coefficient and that coefficient involves padic_.

	Compare1() in sort.c sets AT.SortPadicMode to indicate which term(s) carry a
	padic_ coefficient record:
	- 3: both terms end with a proper padic_ record
	- 1: only *ps1 has padic_, *ps2 has a rational coefficient
	- 2: only *ps2 has padic_, *ps1 has a rational coefficient

	AddWithPadic computes the p-adic sum and rewrites *ps1 to contain exactly one
	padic_ record and a unit coefficient 1/1. It tries to reuse term space in
	place and otherwise allocates in the sort buffer.
*/
int AddWithPadic(PHEAD WORD **ps1, WORD **ps2)
{
	GETBIDENTITY
	SORTING *S = AT.SS;
	PADIC_AUX *aux;
	WORD *coef1, *coef2, size1, size2, *fun1, *fun2, *fun3;
	WORD *s1, *s2, *t1, *t2, i, j, jj;

	if ( !PadicActive ) return(0);
	aux = PadicAux;

	s1 = *ps1;
	s2 = *ps2;
	coef1 = s1 + *s1; size1 = coef1[-1]; coef1 -= ABS(size1);
	coef2 = s2 + *s2; size2 = coef2[-1]; coef2 -= ABS(size2);

	if ( AT.SortPadicMode == 3 ) {
		/* Both coefficients are padic_: unpack and apply external +/- sign. */
		fun1 = s1+1; while ( fun1 < coef1 && fun1[0] != PADICFUN ) fun1 += fun1[1];
		fun2 = s2+1; while ( fun2 < coef2 && fun2[0] != PADICFUN ) fun2 += fun2[1];
		UnpackPadic(aux->p1,fun1);
		if ( size1 < 0 ) padic_neg(aux->p1,aux->p1,PadicContext);
		UnpackPadic(aux->p2,fun2);
		if ( size2 < 0 ) padic_neg(aux->p2,aux->p2,PadicContext);
	}
	else if ( AT.SortPadicMode == 1 ) {
		/* First coefficient is padic_, second is rational. */
		fun1 = s1+1; while ( fun1 < coef1 && fun1[0] != PADICFUN ) fun1 += fun1[1];
		UnpackPadic(aux->p1,fun1);
		if ( size1 < 0 ) padic_neg(aux->p1,aux->p1,PadicContext);
		FormRatToFmpq(aux->q1,(UWORD *)coef2,size2);
		padic_set_fmpq(aux->p2,aux->q1,PadicContext);
	}
	else if ( AT.SortPadicMode == 2 ) {
		/* Second coefficient is padic_, first is rational. */
		fun2 = s2+1; while ( fun2 < coef2 && fun2[0] != PADICFUN ) fun2 += fun2[1];
		UnpackPadic(aux->p2,fun2);
		if ( size2 < 0 ) padic_neg(aux->p2,aux->p2,PadicContext);
		FormRatToFmpq(aux->q1,(UWORD *)coef1,size1);
		padic_set_fmpq(aux->p1,aux->q1,PadicContext);
	}
	else {
		MLOCK(ErrorMessageLock);
		MesPrint("Illegal value %d for AT.SortPadicMode in AddWithPadic.",AT.SortPadicMode);
		MUNLOCK(ErrorMessageLock);
		Terminate(-1);
		return(0);
	}

	padic_add(aux->p3,aux->p1,aux->p2,PadicContext);
	if ( padic_is_zero(aux->p3) ) {
		/* Terms cancel. */
		*ps1 = *ps2 = 0;
		AT.SortPadicMode = 0;
		return(0);
	}

	fun3 = TermMalloc("AddWithPadic");
	PackPadic(fun3,aux->p3);

	if ( AT.SortPadicMode == 3 ) {
		/* Prefer overwriting an existing padic_ record in-place if it fits. */
		if ( fun1[1] == fun3[1] ) {
Over1:
			i = fun3[1]; t1 = fun1; t2 = fun3; NCOPY(t1,t2,i);
			*t1++ = 1; *t1++ = 1; *t1++ = 3;
			*s1 = t1-s1; goto Finished;
		}
		else if ( fun2[1] == fun3[1] ) {
Over2:
			i = fun3[1]; t1 = fun2; t2 = fun3; NCOPY(t1,t2,i);
			*t1++ = 1; *t1++ = 1; *t1++ = 3;
			*s2 = t1-s2; *ps1 = s2; goto Finished;
		}
		else if ( fun1[1] >= fun3[1] ) goto Over1;
		else if ( fun2[1] >= fun3[1] ) goto Over2;
	}
	else if ( AT.SortPadicMode == 1 ) {
		if ( fun1[1] >= fun3[1] ) goto Over1;
		else if ( fun3[1]+3 <= ABS(size2) ) goto Over2;
	}
	else if ( AT.SortPadicMode == 2 ) {
		if ( fun2[1] >= fun3[1] ) goto Over2;
		else if ( fun3[1]+3 <= ABS(size1) ) goto Over1;
	}

	jj = fun1-s1;
	/* Required space: prefix + padic_ record + coefficient words (1,1,3). */
	j = jj+fun3[1]+3;
	if ( (S->sFill + j) >= S->sTop2 ) {
		/* Make space in the sort buffer and refresh pointers. */
		GarbHand();
		s1 = *ps1;
		fun1 = s1+jj;
	}
	t1 = S->sFill;
	for ( i = 0; i < jj; i++ ) *t1++ = s1[i];
	i = fun3[1]; s1 = fun3; NCOPY(t1,s1,i);
	*t1++ = 1; *t1++ = 1; *t1++ = 3;
	*ps1 = S->sFill;
	**ps1 = t1-*ps1;
	S->sFill = t1;

Finished:
	*ps2 = 0;
	TermFree(fun3,"AddWithPadic");
	AT.SortPadicMode = 0;
	if ( **ps1 > AM.MaxTer/((LONG)(sizeof(WORD))) ) {
		MLOCK(ErrorMessageLock);
		MesPrint("Term too complex after p-adic addition in sort. MaxTermSize = %10l",
			AM.MaxTer/sizeof(WORD));
		MUNLOCK(ErrorMessageLock);
		Terminate(-1);
	}
	return(1);
}
/*
  		#] AddWithPadic :
  		#[ MergeWithPadic :

	Variant of AddWithPadic used during patch merging: it computes the p-adic sum
	and rewrites term1 in-place if possible, otherwise it shifts/overwrites
	memory so that *interm1 points to the updated term.
*/
int MergeWithPadic(PHEAD WORD **interm1, WORD **interm2)
{
	GETBIDENTITY
	PADIC_AUX *aux;
	WORD *coef1, *coef2, size1, size2, *fun1, *fun2, *fun3, *tt;
	WORD jj, *t1, *t2, i, *term1 = *interm1, *term2 = *interm2;
	int retval = 0;

	if ( !PadicActive ) return(0);
	aux = PadicAux;

	coef1 = term1+*term1; size1 = coef1[-1]; coef1 -= ABS(size1);
	coef2 = term2+*term2; size2 = coef2[-1]; coef2 -= ABS(size2);
	if ( AT.SortPadicMode == 3 ) {
		fun1 = term1+1; while ( fun1 < coef1 && fun1[0] != PADICFUN ) fun1 += fun1[1];
		fun2 = term2+1; while ( fun2 < coef2 && fun2[0] != PADICFUN ) fun2 += fun2[1];
		UnpackPadic(aux->p1,fun1);
		if ( size1 < 0 ) padic_neg(aux->p1,aux->p1,PadicContext);
		UnpackPadic(aux->p2,fun2);
		if ( size2 < 0 ) padic_neg(aux->p2,aux->p2,PadicContext);
	}
	else if ( AT.SortPadicMode == 1 ) {
		fun1 = term1+1; while ( fun1 < coef1 && fun1[0] != PADICFUN ) fun1 += fun1[1];
		UnpackPadic(aux->p1,fun1);
		if ( size1 < 0 ) padic_neg(aux->p1,aux->p1,PadicContext);
		FormRatToFmpq(aux->q1,(UWORD *)coef2,size2);
		padic_set_fmpq(aux->p2,aux->q1,PadicContext);
	}
	else if ( AT.SortPadicMode == 2 ) {
		fun2 = term2+1; while ( fun2 < coef2 && fun2[0] != PADICFUN ) fun2 += fun2[1];
		FormRatToFmpq(aux->q1,(UWORD *)coef1,size1);
		padic_set_fmpq(aux->p1,aux->q1,PadicContext);
		UnpackPadic(aux->p2,fun2);
		if ( size2 < 0 ) padic_neg(aux->p2,aux->p2,PadicContext);
	}
	else {
		MLOCK(ErrorMessageLock);
		MesPrint("Illegal value %d for AT.SortPadicMode in MergeWithPadic.",AT.SortPadicMode);
		MUNLOCK(ErrorMessageLock);
		Terminate(-1);
		return(0);
	}

	padic_add(aux->p3,aux->p1,aux->p2,PadicContext);
	if ( padic_is_zero(aux->p3) ) {
		AT.SortPadicMode = 0;
		return(0);
	}

	fun3 = TermMalloc("MergeWithPadic");
	PackPadic(fun3,aux->p3);
		if ( AT.SortPadicMode == 3 ) {
			if ( fun1[1] + ABS(size1) == fun3[1] + 3 ) {
OnTopOf1:
				/* The new (padic_ + 1/1) fits exactly on top of the old suffix. */
				t1 = fun3; t2 = fun1;
				for ( i = 0; i < fun3[1]; i++ ) *t2++ = *t1++;
				*t2++ = 1; *t2++ = 1; *t2++ = 3;
				retval = 1;
			}
			else if ( fun1[1] + ABS(size1) > fun3[1] + 3 ) {
Shift1:
				/* There is slack in term1; shift the tail down and rewrite in place. */
				t2 = term1 + *term1; tt = t2;
				*--t2 = 3; *--t2 = 1; *--t2 = 1;
				t1 = fun3 + fun3[1];
				for ( i = 0; i < fun3[1]; i++ ) *--t2 = *--t1;
			t1 = fun1;
			while ( t1 > term1 ) *--t2 = *--t1;
			*t2 = tt-t2; term1 = t2;
			retval = 1;
			}
			else {
				jj = fun3[1]-fun1[1]+3-ABS(size1);
Over1:
				/* term1 needs to grow: move the start pointer back by jj words. */
				t2 = term1-jj; t1 = term1;
				while ( t1 < fun1 ) *t2++ = *t1++;
				term1 -= jj;
				*term1 += jj;
			for ( i = 0; i < fun3[1]; i++ ) *t2++ = fun3[i];
			*t2++ = 1; *t2++ = 1; *t2++ = 3;
			retval = 1;
		}
	}
	else if ( AT.SortPadicMode == 1 ) {
		if ( fun1[1] + ABS(size1) == fun3[1] + 3 ) goto OnTopOf1;
		else if ( fun1[1] + ABS(size1) > fun3[1] + 3 ) goto Shift1;
		else {
			jj = fun3[1]-fun1[1]+3-ABS(size1);
			goto Over1;
		}
	}
	else {
		if ( fun3[1] + 3 == ABS(size1) ) goto OnTopOf1;
		else if ( fun3[1] + 3 < ABS(size1) ) goto Shift1;
		else {
			jj = fun3[1]+3-ABS(size1);
			goto Over1;
		}
	}
	*interm1 = term1;
	TermFree(fun3,"MergeWithPadic");
	AT.SortPadicMode = 0;
	return(retval);
}
/*
 		#] MergeWithPadic :
  	#] Sorting :
*/
