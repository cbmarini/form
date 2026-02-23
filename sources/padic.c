/** @file padic.c
 *
 *  This file contains the p-adic runtime integration.
 *  The implementation follows the same integration points as float.c:
 *  - a dedicated internal function `padic_` used as coefficient carrier,
 *  - statement support (`ToPadic`) in compiler/executor,
 *  - normalization/sorting helper routines for coefficient arithmetic,
 *  - print support with `Format padicprecision`.
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
  	#[ Runtime state :

	FORM keeps a single active p-adic context per run:
	- prime p
	- precision N
	The context is configured by #StartPadic / #EndPadic.
*/
#define PadicRuntimeActive     AM.PadicRuntimeActive
#define PadicContextInitialized AM.PadicContextInitialized
#define PadicPrime             AM.PadicPrime
#define PadicPrecision         AM.PadicPrecision
#define ActivePadicContext     ((padic_ctx_struct *)(AM.PadicContext))

/*
	Thread-local aux storage. The pointer is stored in AT.padic_aux_.
*/
typedef struct PADIC_AUX_ {
	padic_t p1;
	padic_t p2;
	padic_t p3;
	padic_t p4;
	mpq_t q1;
	mpz_t z1;
	mpz_t z2;
} PADIC_AUX;

/*
  	#] Runtime state :
  	#[ Helpers :
 		#[ GetPadicAux :
*/
static PADIC_AUX *GetPadicAux(void)
{
	GETIDENTITY
	return (PADIC_AUX *)(AT.padic_aux_);
}
/*
 		#] GetPadicAux :
 		#[ InitPadicAux :
*/
static void InitPadicAux(PADIC_AUX *aux, LONG prec)
{
	padic_init2(aux->p1, (slong)prec);
	padic_init2(aux->p2, (slong)prec);
	padic_init2(aux->p3, (slong)prec);
	padic_init2(aux->p4, (slong)prec);
	mpq_init(aux->q1);
	mpz_init(aux->z1);
	mpz_init(aux->z2);
}
/*
 		#] InitPadicAux :
 		#[ ClearPadicAux :
*/
static void ClearPadicAux(PADIC_AUX *aux)
{
	padic_clear(aux->p1);
	padic_clear(aux->p2);
	padic_clear(aux->p3);
	padic_clear(aux->p4);
	mpq_clear(aux->q1);
	mpz_clear(aux->z1);
	mpz_clear(aux->z2);
}
/*
 		#] ClearPadicAux :
 		#[ AllocatePadicAuxForAllThreads :
*/
static int AllocatePadicAuxForAllThreads(void)
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
			ClearPadicAux((PADIC_AUX *)(AB[id]->T.padic_aux_));
			M_free(AB[id]->T.padic_aux_,"AB[id]->T.padic_aux_");
			AB[id]->T.padic_aux_ = 0;
		}
		aux = (PADIC_AUX *)Malloc1(sizeof(PADIC_AUX),"AB[id]->T.padic_aux_");
		if ( aux == 0 ) return(-1);
		InitPadicAux(aux,PadicPrecision);
		AB[id]->T.padic_aux_ = (void *)aux;
	}
#else
	PADIC_AUX *aux;
	if ( AT.padic_aux_ ) {
		ClearPadicAux((PADIC_AUX *)(AT.padic_aux_));
		M_free(AT.padic_aux_,"AT.padic_aux_");
		AT.padic_aux_ = 0;
	}
	aux = (PADIC_AUX *)Malloc1(sizeof(PADIC_AUX),"AT.padic_aux_");
	if ( aux == 0 ) return(-1);
	InitPadicAux(aux,PadicPrecision);
	AT.padic_aux_ = (void *)aux;
#endif
	return(0);
}
/*
 		#] AllocatePadicAuxForAllThreads :
 		#[ ClearPadicAuxForAllThreads :
*/
static void ClearPadicAuxForAllThreads(void)
{
#ifdef WITHPTHREADS
	int id, totnum;
	totnum = AM.totalnumberofthreads;
#ifdef WITHSORTBOTS
	totnum = MaX(2 * AM.totalnumberofthreads - 3, AM.totalnumberofthreads);
#endif
	for ( id = 0; id < totnum; id++ ) {
		if ( AB[id]->T.padic_aux_ ) {
			ClearPadicAux((PADIC_AUX *)(AB[id]->T.padic_aux_));
			M_free(AB[id]->T.padic_aux_,"AB[id]->T.padic_aux_");
			AB[id]->T.padic_aux_ = 0;
		}
	}
#else
	if ( AT.padic_aux_ ) {
		ClearPadicAux((PADIC_AUX *)(AT.padic_aux_));
		M_free(AT.padic_aux_,"AT.padic_aux_");
		AT.padic_aux_ = 0;
	}
#endif
}
/*
 		#] ClearPadicAuxForAllThreads :
 		#[ WriteSmallArg :
*/
static WORD *WriteSmallArg(WORD *t, WORD value)
{
	*t++ = -SNUMBER;
	*t++ = value;
	return(t);
}
/*
 		#] WriteSmallArg :
 		#[ WriteLongArg :

	Writes a signed LONG argument in either compact -SNUMBER form or as
	a two-word long integer argument (same internal style used in float.c
	for signed long exponents).
*/
static WORD *WriteLongArg(WORD *t, LONG value)
{
	ULONG x;
	if ( value >= WORD_MIN_VALUE && value <= WORD_MAX_VALUE ) {
		return(WriteSmallArg(t,(WORD)value));
	}
	x = ( value < 0 ) ? (ULONG)(-value) : (ULONG)value;
	*t++ = ARGHEAD + 6;
	*t++ = 0;
	FILLARG(t);
	*t++ = 6;
	*t++ = (UWORD)x;
	*t++ = (UWORD)(x >> BITSINWORD);
	*t++ = 1;
	*t++ = 0;
	*t++ = (value < 0) ? -5 : 5;
	return(t);
}
/*
 		#] WriteLongArg :
 		#[ ReadLongArg :

	Reads a signed LONG argument written by WriteLongArg.
	Returns 0 on failure.
*/
static WORD *ReadLongArg(WORD *t, LONG *value)
{
	if ( *t == -SNUMBER ) {
		*value = (LONG)t[1];
		return(t+2);
	}
	if ( *t == ARGHEAD+6
	 && ABS(t[ARGHEAD+5]) == 5
	 && t[ARGHEAD+3] == 1
	 && t[ARGHEAD+4] == 0 ) {
		ULONG x = ((ULONG)(UWORD)t[ARGHEAD+2] << BITSINWORD) + (UWORD)t[ARGHEAD+1];
		*value = (t[ARGHEAD+5] < 0) ? -(LONG)x : (LONG)x;
		return(t + *t);
	}
	return(0);
}
/*
 		#] ReadLongArg :
 		#[ FormRatToMpq :

	Converts the internal Form coefficient representation (formrat/ratsize)
	to a GMP rational.
*/
static void FormRatToMpq(mpq_t result, UWORD *formrat, WORD ratsize)
{
	WORD nnum, nden;
	UWORD *num, *den;
	int sign = 1;
	mpz_t znum, zden;

	mpz_init(znum);
	mpz_init(zden);

	if ( ratsize < 0 ) {
		ratsize = -ratsize;
		sign = -1;
	}
	nnum = nden = (ratsize-1)/2;
	num = formrat;
	den = formrat + nnum;
	while ( nnum > 0 && num[nnum-1] == 0 ) nnum--;
	while ( nden > 0 && den[nden-1] == 0 ) nden--;

	if ( nnum > 0 ) mpz_import(znum,(size_t)nnum,-1,sizeof(UWORD),0,0,num);
	else            mpz_set_ui(znum,0UL);
	if ( nden > 0 ) mpz_import(zden,(size_t)nden,-1,sizeof(UWORD),0,0,den);
	else            mpz_set_ui(zden,1UL);

	if ( sign < 0 ) mpz_neg(znum,znum);

	mpq_set_num(result,znum);
	mpq_set_den(result,zden);
	mpq_canonicalize(result);

	mpz_clear(zden);
	mpz_clear(znum);
}
/*
 		#] FormRatToMpq :
 		#[ WriteMpzArg :
*/
static WORD *WriteMpzArg(WORD *t, mpz_t z, const char *who)
{
	LONG small;
	int sign;
	mpz_t zabs;
	size_t count = 0;
	UWORD *limbs = 0;
	WORD nnum, i;

	if ( mpz_fits_slong_p(z) ) {
		small = (LONG)mpz_get_si(z);
		if ( small >= WORD_MIN_VALUE && small <= WORD_MAX_VALUE ) {
			return(WriteSmallArg(t,(WORD)small));
		}
	}

	sign = mpz_sgn(z);
	if ( sign == 0 ) return(WriteSmallArg(t,0));

	mpz_init(zabs);
	mpz_set(zabs,z);
	if ( sign < 0 ) mpz_neg(zabs,zabs);

	count = (size_t)((mpz_sizeinbase(zabs,2) + BITSINWORD - 1) / BITSINWORD);
	if ( count > (size_t)WORD_MAX_VALUE ) {
		mpz_clear(zabs);
		MLOCK(ErrorMessageLock);
		MesPrint("p-adic internal overflow in %s.",who);
		MUNLOCK(ErrorMessageLock);
		Terminate(-1);
	}

	limbs = (UWORD *)Malloc1(count*sizeof(UWORD),who);
	if ( limbs == 0 ) {
		mpz_clear(zabs);
		MLOCK(ErrorMessageLock);
		MesPrint("Fatal error in Malloc1 call in %s.",who);
		MUNLOCK(ErrorMessageLock);
		Terminate(-1);
	}

	mpz_export(limbs,&count,-1,sizeof(UWORD),0,0,zabs);
	mpz_clear(zabs);
	if ( count == 0 ) {
		M_free(limbs,who);
		return(WriteSmallArg(t,0));
	}
	if ( count > (size_t)WORD_MAX_VALUE ) {
		M_free(limbs,who);
		MLOCK(ErrorMessageLock);
		MesPrint("p-adic internal overflow in %s.",who);
		MUNLOCK(ErrorMessageLock);
		Terminate(-1);
	}

	nnum = (WORD)count;
	*t++ = ARGHEAD + 2*nnum + 2;
	*t++ = 0;
	FILLARG(t);
	*t++ = 2*nnum + 2;
	for ( i = 0; i < nnum; i++ ) *t++ = (WORD)limbs[i];
	*t++ = 1;
	for ( i = 1; i < nnum; i++ ) *t++ = 0;
	*t++ = ( sign < 0 ) ? -(2*nnum+1) : (2*nnum+1);

	M_free(limbs,who);
	return(t);
}
/*
 		#] WriteMpzArg :
 		#[ ReadMpzArg :
*/
static WORD *ReadMpzArg(WORD *f, WORD *fstop, mpz_t z)
{
	WORD nnum, i, signcode;
	if ( f >= fstop ) return(0);

	if ( *f == -SNUMBER ) {
		mpz_set_si(z,(slong)(f[1]));
		return(f+2);
	}

	if ( *f <= 0 || f + *f > fstop ) return(0);
	signcode = f[*f-1];
	if ( signcode == 0 ) return(0);
	nnum = (WORD)((ABS(signcode)-1)/2);
	if ( nnum <= 0 ) return(0);
	if ( ABS(signcode) != 2*nnum + 1 ) return(0);
	if ( *f != ARGHEAD + 2*nnum + 2 ) return(0);
	if ( f[ARGHEAD] != 2*nnum + 2 ) return(0);
	if ( f[ARGHEAD+1+nnum] != 1 ) return(0);
	for ( i = 1; i < nnum; i++ ) {
		if ( f[ARGHEAD+1+nnum+i] != 0 ) return(0);
	}

	mpz_import(z,(size_t)nnum,-1,sizeof(UWORD),0,0,(UWORD *)(f+ARGHEAD+1));
	if ( signcode < 0 ) mpz_neg(z,z);
	return(f + *f);
}
/*
 		#] ReadMpzArg :
 		#[ ContextMismatchError :
*/
static void ContextMismatchError(LONG N)
{
	MLOCK(ErrorMessageLock);
	MesPrint("Incompatible p-adic precision in padic_ coefficient: found N=%l, active context uses p=%l, N=%l.",
		N,PadicPrime,PadicPrecision);
	MesPrint("Use %#StartPadic with matching parameters before combining these terms.");
	MUNLOCK(ErrorMessageLock);
	Terminate(-1);
}
/*
 		#] ContextMismatchError :
  	#[ Internal p-adic function format :
 		#[ UnpackPadic :

	The internal `padic_` representation stores the FLINT triplet (u,v,N):

	  padic_(v, N, u)

	where u is encoded as a normal Form integer argument (either -SNUMBER
	or a long integer n/1 argument, exactly as in float_ limb packing).
*/
static int UnpackPadic(PADIC_AUX *aux, padic_t out, WORD *fun)
{
	WORD *f, *fstop;
	LONG v, N;

	if ( !PadicRuntimeActive || !PadicContextInitialized ) {
		MLOCK(ErrorMessageLock);
		MesPrint("Illegal attempt at using a padic_ function without proper startup.");
		MesPrint("Please use %#StartPadic <p>,N=<N> first.");
		MUNLOCK(ErrorMessageLock);
		Terminate(-1);
	}

	if ( TestPadic(fun) == 0 ) return(-1);

	f = fun + FUNHEAD;
	fstop = fun + fun[1];

	f = ReadLongArg(f,&v);
	if ( f == 0 ) return(-1);
	f = ReadLongArg(f,&N);
	if ( f == 0 ) return(-1);
	if ( N != PadicPrecision ) {
		ContextMismatchError(N);
		return(-1);
	}

	f = ReadMpzArg(f,fstop,aux->z1);
	if ( f == 0 ) return(-1);
	if ( f != fstop ) return(-1);
	fmpz_set_mpz(padic_unit(out),aux->z1);
	padic_val(out) = (slong)v;
	padic_prec(out) = (slong)N;
	padic_reduce(out,ActivePadicContext);
	return(0);
}
/*
 		#] UnpackPadic :
 		#[ PackPadic :
*/
static int PackPadic(PADIC_AUX *aux, WORD *fun, padic_t in)
{
	WORD *t;
	LONG v, N;

	padic_set(aux->p4,in,ActivePadicContext);
	padic_reduce(aux->p4,ActivePadicContext);
	v = (LONG)padic_val(aux->p4);
	N = (LONG)padic_prec(aux->p4);

	fmpz_get_mpz(aux->z1,padic_unit(aux->p4));

	t = fun;
	*t++ = PADICFUN;
	t++;
	FILLFUN(t);

	t = WriteLongArg(t,v);
	t = WriteLongArg(t,N);
	t = WriteMpzArg(t,aux->z1,"PackPadic(unit)");
	fun[1] = t - fun;
	return(fun[1]);
}
/*
 		#] PackPadic :
 		#] Internal p-adic function format :
  	#] Helpers :
  	#[ Runtime lifecycle :
 		#[ PadicIsActive :
*/
int PadicIsActive(void)
{
	return(PadicRuntimeActive);
}
/*
 		#] PadicIsActive :
 		#[ PadicIsPrime :
*/
int PadicIsPrime(LONG p)
{
	int prime;
	fmpz_t z;
	if ( p <= 1 ) return(0);
	fmpz_init(z);
	fmpz_set_si(z,(slong)p);
	prime = fmpz_is_prime(z);
	fmpz_clear(z);
	return(prime == 1);
}
/*
 		#] PadicIsPrime :
 		#[ StartPadicSystem :
*/
int StartPadicSystem(LONG p, LONG N)
{
	fmpz_t prime;
	if ( p <= 1 || N <= 0 ) return(1);
	if ( PadicRuntimeActive ) {
		ClearPadicSystem();
	}
	if ( AM.PadicContext == 0 ) {
		AM.PadicContext = Malloc1(sizeof(padic_ctx_struct),"PadicContext");
		if ( AM.PadicContext == 0 ) return(1);
	}
	if ( PadicContextInitialized ) {
		padic_ctx_clear(ActivePadicContext);
		PadicContextInitialized = 0;
	}
	PadicPrime = p;
	PadicPrecision = N;
	fmpz_init(prime);
	fmpz_set_si(prime,(slong)p);
	padic_ctx_init(ActivePadicContext,prime,0,(slong)N,PADIC_SERIES);
	fmpz_clear(prime);
	PadicContextInitialized = 1;
	PadicRuntimeActive = 1;

	if ( AllocatePadicAuxForAllThreads() ) {
		MLOCK(ErrorMessageLock);
		MesPrint("Failed to initialize p-adic thread-local buffers.");
		MUNLOCK(ErrorMessageLock);
		ClearPadicSystem();
		return(1);
	}
	return(0);
}
/*
 		#] StartPadicSystem :
 		#[ ClearPadicSystem :
*/
void ClearPadicSystem(void)
{
	ClearPadicAuxForAllThreads();
	if ( PadicContextInitialized ) {
		padic_ctx_clear(ActivePadicContext);
		PadicContextInitialized = 0;
	}
	if ( AM.PadicContext ) {
		M_free(AM.PadicContext,"PadicContext");
		AM.PadicContext = 0;
	}
	PadicRuntimeActive = 0;
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
  	#] Runtime lifecycle :
  	#[ Validation and conversion :
 		#[ TestPadic :
*/
int TestPadic(WORD *fun)
{
	WORD *f, *fstop;
	WORD nargs;
	LONG v, N;
	mpz_t z;

	f = fun + FUNHEAD;
	fstop = fun + fun[1];
	nargs = 0;
	while ( f < fstop ) {
		nargs++;
		NEXTARG(f);
	}
	if ( nargs != 3 ) return(0);

	f = fun + FUNHEAD;
	f = ReadLongArg(f,&v);
	if ( f == 0 ) return(0);
	f = ReadLongArg(f,&N);
	if ( f == 0 || N <= 0 ) return(0);
	mpz_init(z);
	f = ReadMpzArg(f,fstop,z);
	mpz_clear(z);
	if ( f == 0 || f != fstop ) return(0);
	return(1);
}
/*
 		#] TestPadic :
 		#[ RatToPadicFun :
*/
int RatToPadicFun(PHEAD WORD *outfun, UWORD *formrat, WORD nrat)
{
	PADIC_AUX *aux;
	if ( !PadicRuntimeActive ) return(-1);
	aux = (PADIC_AUX *)(AT.padic_aux_);
	if ( aux == 0 ) return(-1);
	FormRatToMpq(aux->q1,formrat,nrat);
	padic_set_mpq(aux->p1,aux->q1,ActivePadicContext);
	PackPadic(aux,outfun,aux->p1);
	return(0);
}
/*
 		#] RatToPadicFun :
 		#[ MulRatToPadic :
*/
int MulRatToPadic(PHEAD WORD *outfun, WORD *infun, UWORD *formrat, WORD nrat)
{
	PADIC_AUX *aux;
	if ( !PadicRuntimeActive ) return(-1);
	aux = (PADIC_AUX *)(AT.padic_aux_);
	if ( aux == 0 ) return(-1);
	if ( UnpackPadic(aux,aux->p1,infun) ) return(-1);
	FormRatToMpq(aux->q1,formrat,nrat);
	padic_set_mpq(aux->p2,aux->q1,ActivePadicContext);
	padic_mul(aux->p3,aux->p1,aux->p2,ActivePadicContext);
	PackPadic(aux,outfun,aux->p3);
	return(0);
}
/*
 		#] MulRatToPadic :
 		#[ MulPadics :
*/
int MulPadics(PHEAD WORD *fun3, WORD *fun1, WORD *fun2)
{
	PADIC_AUX *aux;
	if ( !PadicRuntimeActive ) return(-1);
	aux = (PADIC_AUX *)(AT.padic_aux_);
	if ( aux == 0 ) return(-1);
	if ( UnpackPadic(aux,aux->p1,fun1) ) return(-1);
	if ( UnpackPadic(aux,aux->p2,fun2) ) return(-1);
	padic_mul(aux->p3,aux->p1,aux->p2,ActivePadicContext);
	PackPadic(aux,fun3,aux->p3);
	return(0);
}
/*
 		#] MulPadics :
 		#[ DivPadics :
*/
int DivPadics(PHEAD WORD *fun3, WORD *fun1, WORD *fun2)
{
	PADIC_AUX *aux;
	if ( !PadicRuntimeActive ) return(-1);
	aux = (PADIC_AUX *)(AT.padic_aux_);
	if ( aux == 0 ) return(-1);
	if ( UnpackPadic(aux,aux->p1,fun1) ) return(-1);
	if ( UnpackPadic(aux,aux->p2,fun2) ) return(-1);
	if ( padic_is_zero(aux->p2) ) {
		MLOCK(ErrorMessageLock);
		MesPrint("Division by zero in p-adic arithmetic.");
		MUNLOCK(ErrorMessageLock);
		Terminate(-1);
		return(-1);
	}
	padic_div(aux->p3,aux->p1,aux->p2,ActivePadicContext);
	PackPadic(aux,fun3,aux->p3);
	return(0);
}
/*
 		#] DivPadics :
 		#[ MpqToFormRat :
*/
static int MpqToFormRat(UWORD *out, WORD *nratout, mpq_t q)
{
	int sign;
	size_t nnum, nden, n, i, count;
	mpz_t znum;

	sign = mpq_sgn(q);
	if ( sign == 0 ) {
		*nratout = 3;
		if ( out != 0 ) {
			out[0] = 0;
			out[1] = 1;
			out[2] = 3;
		}
		return(0);
	}

	nnum = (mpz_sgn(mpq_numref(q)) == 0) ? 0 :
		(size_t)((mpz_sizeinbase(mpq_numref(q),2) + BITSINWORD - 1) / BITSINWORD);
	nden = (mpz_sgn(mpq_denref(q)) == 0) ? 0 :
		(size_t)((mpz_sizeinbase(mpq_denref(q),2) + BITSINWORD - 1) / BITSINWORD);
	if ( nden == 0 ) return(-1);
	n = ( nnum > nden ) ? nnum : nden;
	if ( n > (size_t)((WORD_MAX_VALUE-1)/2) ) return(-1);

	*nratout = (WORD)(2*n + 1);
	if ( sign < 0 ) *nratout = -*nratout;
	if ( out == 0 ) return(0);

	for ( i = 0; i < 2*n; i++ ) out[i] = 0;

	mpz_init(znum);
	mpz_set(znum,mpq_numref(q));
	if ( mpz_sgn(znum) < 0 ) mpz_neg(znum,znum);
	if ( nnum > 0 ) {
		count = nnum;
		mpz_export(out,&count,-1,sizeof(UWORD),0,0,znum);
	}
	if ( nden > 0 ) {
		count = nden;
		mpz_export(out+n,&count,-1,sizeof(UWORD),0,0,mpq_denref(q));
	}
	mpz_clear(znum);
	out[2*n] = (UWORD)ABS(*nratout);
	return(0);
}
/*
 		#] MpqToFormRat :
 		#[ PadicReconstructToMpq :

	Reconstructs a small rational from a reduced p-adic value using FLINT's
	rational reconstruction and applies the p-adic valuation afterwards.
*/
static int PadicReconstructToMpq(mpq_t out, padic_t in)
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

	fmpz_pow_ui(modulus,ActivePadicContext->p,(ulong)modexp);
	fmpz_mod(residue,padic_unit(in),modulus);
	ok = fmpq_reconstruct_fmpz(recon,residue,modulus);
	if ( !ok ) goto ClearAndReturn;

	if ( v > 0 ) {
		fmpz_pow_ui(ppower,ActivePadicContext->p,(ulong)v);
		fmpq_mul_fmpz(recon,recon,ppower);
	}
	else if ( v < 0 ) {
		fmpz_pow_ui(ppower,ActivePadicContext->p,(ulong)(-v));
		fmpq_div_fmpz(recon,recon,ppower);
	}

	fmpq_get_mpq(out,recon);

ClearAndReturn:
	fmpq_clear(recon);
	fmpz_clear(ppower);
	fmpz_clear(modulus);
	fmpz_clear(residue);
	return(ok ? 0 : -1);
}
/*
 		#] PadicReconstructToMpq :
  	#] Validation and conversion :
  	#[ Printing :
 		#[ PrintPadic :
*/
int PrintPadic(WORD *fun,int numdigits)
{
	PADIC_AUX *aux;
	char *flint_string;
	size_t n;
	int digits = (int)PadicPrecision;

	if ( !PadicRuntimeActive || !PadicContextInitialized ) return(0);
	aux = GetPadicAux();
	if ( aux == 0 ) return(0);
	if ( UnpackPadic(aux,aux->p1,fun) ) return(0);

	if ( numdigits > 0 && numdigits < digits ) digits = numdigits;

	if ( digits == (int)PadicPrecision ) {
		flint_string = padic_get_str(0,aux->p1,ActivePadicContext);
	}
	else {
		padic_ctx_t short_ctx;
		padic_t short_x;
		padic_ctx_init(short_ctx,ActivePadicContext->p,0,(slong)digits,PADIC_SERIES);
		padic_init2(short_x,digits);
		padic_get_mpq(aux->q1,aux->p1,ActivePadicContext);
		padic_set_mpq(short_x,aux->q1,short_ctx);
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
		if ( AO.padicspace == 0 || AO.padicsize <= (LONG)new_len ) {
			if ( AO.padicspace ) M_free(AO.padicspace,"padicspace");
			AO.padicsize = (LONG)new_len + 32;
			AO.padicspace = (UBYTE *)Malloc1((size_t)AO.padicsize,"padicspace");
		}
		((char *)AO.padicspace)[0] = '(';
		memcpy((char *)AO.padicspace + 1, flint_string, n);
		memcpy((char *)AO.padicspace + 1 + n, tail, tail_len);
		((char *)AO.padicspace)[1 + n + tail_len] = ')';
		((char *)AO.padicspace)[new_len] = 0;
		flint_free(flint_string);
		return((int)new_len);
	}

	n += 2;
	if ( AO.padicspace == 0 || AO.padicsize <= (LONG)n ) {
		if ( AO.padicspace ) M_free(AO.padicspace,"padicspace");
		AO.padicsize = (LONG)n + 32;
		AO.padicspace = (UBYTE *)Malloc1((size_t)AO.padicsize,"padicspace");
	}
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
*/
int CoToPadic(UBYTE *s)
{
	if ( !PadicRuntimeActive ) {
		MesPrint("&Illegal attempt to convert to padic_ without activating p-adic numbers.");
		MesPrint("&Forgotten %#startpadic instruction?");
		return(1);
	}
	while ( *s == ' ' || *s == ',' || *s == '\t' ) s++;
	if ( *s ) {
		MesPrint("&Illegal argument(s) in Topadic statement: '%s'",s);
		return(1);
	}
	Add2Com(TYPETOPADIC);
	return(0);
}
/*
 		#] CoToPadic :
 		#[ CoPadicToRat :
*/
int CoPadicToRat(UBYTE *s)
{
	if ( !PadicRuntimeActive ) {
		MesPrint("&Illegal attempt to convert from padic_ without activating p-adic numbers.");
		MesPrint("&Forgotten %#startpadic instruction?");
		return(1);
	}
	while ( *s == ' ' || *s == ',' || *s == '\t' ) s++;
	if ( *s ) {
		MesPrint("&Illegal argument(s) in PadicToRat statement: '%s'",s);
		return(1);
	}
	Add2Com(TYPETOPADICTORAT);
	return(0);
}
/*
 		#] CoPadicToRat :
 		#[ ToPadic :
*/
int ToPadic(PHEAD WORD *term, WORD level)
{
	GETBIDENTITY
	PADIC_AUX *aux;
	WORD *t, *scan, *tstop, nsize, ncoef;

	if ( !PadicRuntimeActive ) return(1);
	aux = GetPadicAux();
	if ( aux == 0 ) return(1);

	t = term + *term;
	ncoef = t[-1];
	nsize = ABS(ncoef);
	tstop = t - nsize;

	if ( nsize == 3 && tstop[0] == 1 && tstop[1] == 1 ) {
		scan = term + 1;
		while ( scan < tstop ) {
			if ( *scan == PADICFUN && scan + scan[1] == tstop && TestPadic(scan) ) {
				return(Generator(BHEAD term,level));
			}
			scan += scan[1];
		}
	}

	FormRatToMpq(aux->q1,(UWORD *)tstop,ncoef);
	padic_set_mpq(aux->p1,aux->q1,ActivePadicContext);
	PackPadic(aux,tstop,aux->p1);
	tstop += tstop[1];
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
*/
int PadicToRat(PHEAD WORD *term, WORD level)
{
	GETBIDENTITY
	PADIC_AUX *aux;
	WORD *tstop, *t, *stop, nsize, nsign, ncoef;

	if ( !PadicRuntimeActive ) return(1);
	aux = GetPadicAux();
	if ( aux == 0 ) return(1);

	tstop = term + *term;
	nsize = ABS(tstop[-1]);
	nsign = tstop[-1] < 0 ? -1 : 1;
	tstop -= nsize;
	t = term + 1;
	while ( t < tstop ) {
		if ( *t == PADICFUN && t + t[1] == tstop && TestPadic(t)
		 && nsize == 3 && tstop[0] == 1 && tstop[1] == 1 ) break;
		t += t[1];
	}
	if ( t < tstop ) {
		if ( UnpackPadic(aux,aux->p1,t) ) return(1);
		if ( padic_is_zero(aux->p1) ) return(0);
		if ( PadicReconstructToMpq(aux->q1,aux->p1) ) {
			/*
				No unique small reconstruction exists for the current precision.
				Fall back to FLINT's canonical lift.
			*/
			padic_get_mpq(aux->q1,aux->p1,ActivePadicContext);
		}
		if ( MpqToFormRat(0,&ncoef,aux->q1) ) goto RatFailure;
		stop = (WORD *)(((UBYTE *)term) + AM.MaxTer);
		if ( t + ABS(ncoef) > stop ) {
			MLOCK(ErrorMessageLock);
			MesPrint("Term too complex after p-adic to rational conversion. MaxTermSize = %10l",
				AM.MaxTer/sizeof(WORD));
			MUNLOCK(ErrorMessageLock);
			Terminate(-1);
			return(1);
		}
		if ( MpqToFormRat((UWORD *)t,&ncoef,aux->q1) ) goto RatFailure;
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
*/
int AddWithPadic(PHEAD WORD **ps1, WORD **ps2)
{
	GETBIDENTITY
	SORTING *S = AT.SS;
	PADIC_AUX *aux;
	WORD *coef1, *coef2, size1, size2, *fun1, *fun2, *fun3;
	WORD *s1, *s2, *t1, *t2, i, j, jj;

	aux = GetPadicAux();
	if ( aux == 0 ) return(0);

	s1 = *ps1;
	s2 = *ps2;
	coef1 = s1 + *s1; size1 = coef1[-1]; coef1 -= ABS(size1);
	coef2 = s2 + *s2; size2 = coef2[-1]; coef2 -= ABS(size2);

	if ( AT.SortPadicMode == 3 ) {
		fun1 = s1+1; while ( fun1 < coef1 && fun1[0] != PADICFUN ) fun1 += fun1[1];
		fun2 = s2+1; while ( fun2 < coef2 && fun2[0] != PADICFUN ) fun2 += fun2[1];
		UnpackPadic(aux,aux->p1,fun1);
		if ( size1 < 0 ) padic_neg(aux->p1,aux->p1,ActivePadicContext);
		UnpackPadic(aux,aux->p2,fun2);
		if ( size2 < 0 ) padic_neg(aux->p2,aux->p2,ActivePadicContext);
	}
	else if ( AT.SortPadicMode == 1 ) {
		fun1 = s1+1; while ( fun1 < coef1 && fun1[0] != PADICFUN ) fun1 += fun1[1];
		UnpackPadic(aux,aux->p1,fun1);
		if ( size1 < 0 ) padic_neg(aux->p1,aux->p1,ActivePadicContext);
		FormRatToMpq(aux->q1,(UWORD *)coef2,size2);
		padic_set_mpq(aux->p2,aux->q1,ActivePadicContext);
	}
	else if ( AT.SortPadicMode == 2 ) {
		fun2 = s2+1; while ( fun2 < coef2 && fun2[0] != PADICFUN ) fun2 += fun2[1];
		UnpackPadic(aux,aux->p2,fun2);
		if ( size2 < 0 ) padic_neg(aux->p2,aux->p2,ActivePadicContext);
		FormRatToMpq(aux->q1,(UWORD *)coef1,size1);
		padic_set_mpq(aux->p1,aux->q1,ActivePadicContext);
	}
	else {
		MLOCK(ErrorMessageLock);
		MesPrint("Illegal value %d for AT.SortPadicMode in AddWithPadic.",AT.SortPadicMode);
		MUNLOCK(ErrorMessageLock);
		Terminate(-1);
		return(0);
	}

	padic_add(aux->p3,aux->p1,aux->p2,ActivePadicContext);
	if ( padic_is_zero(aux->p3) ) {
		*ps1 = *ps2 = 0;
		AT.SortPadicMode = 0;
		return(0);
	}

	fun3 = TermMalloc("AddWithPadic");
	PackPadic(aux,fun3,aux->p3);

	if ( AT.SortPadicMode == 3 ) {
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
	j = jj+fun3[1]+3;
	if ( (S->sFill + j) >= S->sTop2 ) {
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
*/
int MergeWithPadic(PHEAD WORD **interm1, WORD **interm2)
{
	GETBIDENTITY
	PADIC_AUX *aux;
	WORD *coef1, *coef2, size1, size2, *fun1, *fun2, *fun3, *tt;
	WORD jj, *t1, *t2, i, *term1 = *interm1, *term2 = *interm2;
	int retval = 0;

	aux = GetPadicAux();
	if ( aux == 0 ) return(0);

	coef1 = term1+*term1; size1 = coef1[-1]; coef1 -= ABS(size1);
	coef2 = term2+*term2; size2 = coef2[-1]; coef2 -= ABS(size2);
	if ( AT.SortPadicMode == 3 ) {
		fun1 = term1+1; while ( fun1 < coef1 && fun1[0] != PADICFUN ) fun1 += fun1[1];
		fun2 = term2+1; while ( fun2 < coef2 && fun2[0] != PADICFUN ) fun2 += fun2[1];
		UnpackPadic(aux,aux->p1,fun1);
		if ( size1 < 0 ) padic_neg(aux->p1,aux->p1,ActivePadicContext);
		UnpackPadic(aux,aux->p2,fun2);
		if ( size2 < 0 ) padic_neg(aux->p2,aux->p2,ActivePadicContext);
	}
	else if ( AT.SortPadicMode == 1 ) {
		fun1 = term1+1; while ( fun1 < coef1 && fun1[0] != PADICFUN ) fun1 += fun1[1];
		UnpackPadic(aux,aux->p1,fun1);
		if ( size1 < 0 ) padic_neg(aux->p1,aux->p1,ActivePadicContext);
		FormRatToMpq(aux->q1,(UWORD *)coef2,size2);
		padic_set_mpq(aux->p2,aux->q1,ActivePadicContext);
	}
	else if ( AT.SortPadicMode == 2 ) {
		fun2 = term2+1; while ( fun2 < coef2 && fun2[0] != PADICFUN ) fun2 += fun2[1];
		FormRatToMpq(aux->q1,(UWORD *)coef1,size1);
		padic_set_mpq(aux->p1,aux->q1,ActivePadicContext);
		UnpackPadic(aux,aux->p2,fun2);
		if ( size2 < 0 ) padic_neg(aux->p2,aux->p2,ActivePadicContext);
	}
	else {
		MLOCK(ErrorMessageLock);
		MesPrint("Illegal value %d for AT.SortPadicMode in MergeWithPadic.",AT.SortPadicMode);
		MUNLOCK(ErrorMessageLock);
		Terminate(-1);
		return(0);
	}

	padic_add(aux->p3,aux->p1,aux->p2,ActivePadicContext);
	if ( padic_is_zero(aux->p3) ) {
		AT.SortPadicMode = 0;
		return(0);
	}

	fun3 = TermMalloc("MergeWithPadic");
	PackPadic(aux,fun3,aux->p3);
	if ( AT.SortPadicMode == 3 ) {
		if ( fun1[1] + ABS(size1) == fun3[1] + 3 ) {
OnTopOf1:
			t1 = fun3; t2 = fun1;
			for ( i = 0; i < fun3[1]; i++ ) *t2++ = *t1++;
			*t2++ = 1; *t2++ = 1; *t2++ = 3;
			retval = 1;
		}
		else if ( fun1[1] + ABS(size1) > fun3[1] + 3 ) {
Shift1:
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
*/

/*
  	#] Sorting :
*/
