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
#include <flint/padic.h>

/*
	#[ Runtime state :

	FORM keeps a single active p-adic context per run:
	- prime p
	- precision N
	The context is configured by #StartPadic / #EndPadic.
*/

static int PadicRuntimeActive = 0;
static int PadicPrimeInitialized = 0;
static int PadicContextInitialized = 0;
static LONG PadicPrime = 0;
static LONG PadicPrecision = 0;

static fmpz_t PadicPrimeFmpz;
static padic_ctx_t PadicContext;

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
*/

static PADIC_AUX *GetPadicAux(void)
{
	GETIDENTITY
	return (PADIC_AUX *)(AT.padic_aux_);
}

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

static WORD *WriteSmallArg(WORD *t, WORD value)
{
	*t++ = -SNUMBER;
	*t++ = value;
	return(t);
}

/*
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
	Converts packed limb arguments of the form (-SNUMBER, value) repeated
	`nlimbs` times into an mpz. `limbs` points to the first -SNUMBER tag.
*/
static void PadicLimbsToMpz(mpz_t z, WORD *limbs, WORD nlimbs)
{
	WORD i;
	mpz_set_ui(z,0UL);
	for ( i = nlimbs-1; i >= 0; i-- ) {
		mpz_mul_2exp(z,z,BITSINWORD);
		mpz_add_ui(z,z,(UWORD)limbs[2*i+1]);
	}
}

/*
	Exports a non-negative mpz into little-endian UWORD limbs.
*/
static UWORD *MpzToPadicLimbs(mpz_t z, WORD *nlimbs, const char *who)
{
	size_t count = 0;
	UWORD *out = 0;
	if ( mpz_sgn(z) == 0 ) {
		*nlimbs = 0;
		return(0);
	}
	count = (size_t)((mpz_sizeinbase(z,2) + BITSINWORD - 1) / BITSINWORD);
	if ( count > (size_t)WORD_MAX_VALUE ) {
		MLOCK(ErrorMessageLock);
		MesPrint("p-adic internal overflow in %s.",who);
		MUNLOCK(ErrorMessageLock);
		Terminate(-1);
	}
	out = (UWORD *)Malloc1(count*sizeof(UWORD),who);
	if ( out == 0 ) {
		MLOCK(ErrorMessageLock);
		MesPrint("Fatal error in Malloc1 call in %s.",who);
		MUNLOCK(ErrorMessageLock);
		Terminate(-1);
	}
	mpz_export(out,&count,-1,sizeof(UWORD),0,0,z);
	*nlimbs = (WORD)count;
	return(out);
}

static void ContextMismatchError(LONG p, LONG N)
{
	MLOCK(ErrorMessageLock);
	MesPrint("Incompatible p-adic context in padic_ coefficient: found p=%l, N=%l, active context uses p=%l, N=%l.",
		p,N,PadicPrime,PadicPrecision);
	MesPrint("Use %#StartPadic with matching parameters before combining these terms.");
	MUNLOCK(ErrorMessageLock);
	Terminate(-1);
}

/*
	#[ Internal p-adic function format :

	The internal `padic_` representation stores a context tag and a rational:

	  padic_(p, N, sign, nnum, num_limb_1, ..., num_limb_nnum,
	             nden, den_limb_1, ..., den_limb_nden)

	where the limbs are raw WORD values interpreted as unsigned base 2^BITSINWORD
	digits in little-endian order.
*/

static int UnpackPadic(PADIC_AUX *aux, padic_t out, WORD *fun)
{
	WORD *f, *fstop;
	LONG p, N;
	WORD sign, nnum, nden;

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

	f = ReadLongArg(f,&p);
	if ( f == 0 ) return(-1);
	f = ReadLongArg(f,&N);
	if ( f == 0 ) return(-1);
	if ( p != PadicPrime || N != PadicPrecision ) {
		ContextMismatchError(p,N);
		return(-1);
	}

	sign = f[1];
	f += 2;
	nnum = f[1];
	f += 2;

	PadicLimbsToMpz(aux->z1,f,nnum);
	if ( sign < 0 ) mpz_neg(aux->z1,aux->z1);
	f += 2*nnum;

	nden = f[1];
	f += 2;
	PadicLimbsToMpz(aux->z2,f,nden);
	f += 2*nden;

	if ( f != fstop ) return(-1);
	if ( mpz_sgn(aux->z2) <= 0 ) return(-1);

	mpq_set_num(aux->q1,aux->z1);
	mpq_set_den(aux->q1,aux->z2);
	mpq_canonicalize(aux->q1);
	padic_set_mpq(out,aux->q1,PadicContext);
	return(0);
}

static int PackPadic(PADIC_AUX *aux, WORD *fun, padic_t in)
{
	WORD *t;
	WORD sign, nnum, nden, i;
	UWORD *numlimbs = 0, *denlimbs = 0;

	padic_get_mpq(aux->q1,in,PadicContext);
	mpz_set(aux->z1,mpq_numref(aux->q1));
	mpz_set(aux->z2,mpq_denref(aux->q1));

	sign = (WORD)mpz_sgn(aux->z1);
	if ( sign < 0 ) mpz_neg(aux->z1,aux->z1);

	numlimbs = MpzToPadicLimbs(aux->z1,&nnum,"PackPadic(num)");
	denlimbs = MpzToPadicLimbs(aux->z2,&nden,"PackPadic(den)");
	if ( nden <= 0 ) {
		if ( denlimbs ) M_free(denlimbs,"PackPadic(den)");
		denlimbs = (UWORD *)Malloc1(sizeof(UWORD),"PackPadic(den1)");
		denlimbs[0] = 1;
		nden = 1;
	}
	if ( sign == 0 ) nnum = 0;

	t = fun;
	*t++ = PADICFUN;
	t++;
	FILLFUN(t);

	t = WriteLongArg(t,PadicPrime);
	t = WriteLongArg(t,PadicPrecision);
	t = WriteSmallArg(t,sign);
	t = WriteSmallArg(t,nnum);
	for ( i = 0; i < nnum; i++ ) {
		t = WriteSmallArg(t,(WORD)numlimbs[i]);
	}
	t = WriteSmallArg(t,nden);
	for ( i = 0; i < nden; i++ ) {
		t = WriteSmallArg(t,(WORD)denlimbs[i]);
	}
	fun[1] = t - fun;

	if ( numlimbs ) M_free(numlimbs,"PackPadic(num)");
	if ( denlimbs ) M_free(denlimbs,"PackPadic(den)");
	return(fun[1]);
}

/*
	#] Helpers :
	#[ Runtime lifecycle :
*/

int PadicIsActive(void)
{
	return(PadicRuntimeActive);
}

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

int StartPadicSystem(LONG p, LONG N)
{
	if ( p <= 1 || N <= 0 ) return(1);
	if ( PadicRuntimeActive ) {
		ClearPadicSystem();
	}
	if ( !PadicPrimeInitialized ) {
		fmpz_init(PadicPrimeFmpz);
		PadicPrimeInitialized = 1;
	}
	if ( PadicContextInitialized ) {
		padic_ctx_clear(PadicContext);
		PadicContextInitialized = 0;
	}
	PadicPrime = p;
	PadicPrecision = N;
	fmpz_set_si(PadicPrimeFmpz,(slong)p);
	padic_ctx_init(PadicContext,PadicPrimeFmpz,0,(slong)N,PADIC_SERIES);
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

void ClearPadicSystem(void)
{
	ClearPadicAuxForAllThreads();
	if ( PadicContextInitialized ) {
		padic_ctx_clear(PadicContext);
		PadicContextInitialized = 0;
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
	#] Runtime lifecycle :
	#[ Validation and conversion :
*/

int TestPadic(WORD *fun)
{
	WORD *f, *fstop;
	LONG p, N;
	WORD sign, nnum, nden, i;

	f = fun + FUNHEAD;
	fstop = fun + fun[1];

	f = ReadLongArg(f,&p);
	if ( f == 0 || p <= 1 ) return(0);
	f = ReadLongArg(f,&N);
	if ( f == 0 || N <= 0 ) return(0);

	if ( f >= fstop || *f != -SNUMBER ) return(0);
	sign = f[1];
	if ( sign < -1 || sign > 1 ) return(0);
	f += 2;

	if ( f >= fstop || *f != -SNUMBER ) return(0);
	nnum = f[1];
	if ( nnum < 0 ) return(0);
	f += 2;
	for ( i = 0; i < nnum; i++ ) {
		if ( f >= fstop || *f != -SNUMBER ) return(0);
		f += 2;
	}

	if ( f >= fstop || *f != -SNUMBER ) return(0);
	nden = f[1];
	if ( nden <= 0 ) return(0);
	f += 2;
	for ( i = 0; i < nden; i++ ) {
		if ( f >= fstop || *f != -SNUMBER ) return(0);
		f += 2;
	}
	if ( f != fstop ) return(0);
	if ( sign == 0 && nnum != 0 ) return(0);
	return(1);
}

int RatToPadicFun(PHEAD WORD *outfun, UWORD *formrat, WORD nrat)
{
	PADIC_AUX *aux;
	if ( !PadicRuntimeActive ) return(-1);
	aux = (PADIC_AUX *)(AT.padic_aux_);
	if ( aux == 0 ) return(-1);
	FormRatToMpq(aux->q1,formrat,nrat);
	padic_set_mpq(aux->p1,aux->q1,PadicContext);
	PackPadic(aux,outfun,aux->p1);
	return(0);
}

int MulRatToPadic(PHEAD WORD *outfun, WORD *infun, UWORD *formrat, WORD nrat)
{
	PADIC_AUX *aux;
	if ( !PadicRuntimeActive ) return(-1);
	aux = (PADIC_AUX *)(AT.padic_aux_);
	if ( aux == 0 ) return(-1);
	if ( UnpackPadic(aux,aux->p1,infun) ) return(-1);
	FormRatToMpq(aux->q1,formrat,nrat);
	padic_set_mpq(aux->p2,aux->q1,PadicContext);
	padic_mul(aux->p3,aux->p1,aux->p2,PadicContext);
	PackPadic(aux,outfun,aux->p3);
	return(0);
}

int MulPadics(PHEAD WORD *fun3, WORD *fun1, WORD *fun2)
{
	PADIC_AUX *aux;
	if ( !PadicRuntimeActive ) return(-1);
	aux = (PADIC_AUX *)(AT.padic_aux_);
	if ( aux == 0 ) return(-1);
	if ( UnpackPadic(aux,aux->p1,fun1) ) return(-1);
	if ( UnpackPadic(aux,aux->p2,fun2) ) return(-1);
	padic_mul(aux->p3,aux->p1,aux->p2,PadicContext);
	PackPadic(aux,fun3,aux->p3);
	return(0);
}

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
	padic_div(aux->p3,aux->p1,aux->p2,PadicContext);
	PackPadic(aux,fun3,aux->p3);
	return(0);
}

/*
	#] Validation and conversion :
	#[ Printing :
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
		flint_string = padic_get_str(0,aux->p1,PadicContext);
	}
	else {
		padic_ctx_t short_ctx;
		padic_t short_x;
		padic_ctx_init(short_ctx,PadicPrimeFmpz,0,(slong)digits,PADIC_SERIES);
		padic_init2(short_x,digits);
		padic_get_mpq(aux->q1,aux->p1,PadicContext);
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
		snprintf(tail,sizeof(tail)," + O(%ld^%d)",(long)PadicPrime,digits);
		tail_len = strlen(tail);
		if ( AO.padicspace == 0 || AO.padicsize <= (LONG)(n + tail_len) ) {
			if ( AO.padicspace ) M_free(AO.padicspace,"padicspace");
			AO.padicsize = (LONG)(n + tail_len) + 32;
			AO.padicspace = (UBYTE *)Malloc1((size_t)AO.padicsize,"padicspace");
		}
		strcpy((char *)AO.padicspace,flint_string);
		strcat((char *)AO.padicspace,tail);
		n += tail_len;
		flint_free(flint_string);
		return((int)n);
	}

	if ( AO.padicspace == 0 || AO.padicsize <= (LONG)n ) {
		if ( AO.padicspace ) M_free(AO.padicspace,"padicspace");
		AO.padicsize = (LONG)n + 32;
		AO.padicspace = (UBYTE *)Malloc1((size_t)AO.padicsize,"padicspace");
	}
	strcpy((char *)AO.padicspace,flint_string);
	flint_free(flint_string);
	return((int)n);
}

/*
	#] Printing :
	#[ Compiler/runtime statements :
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

int ToPadic(PHEAD WORD *term, WORD level)
{
	GETBIDENTITY
	PADIC_AUX *aux;
	WORD *t, *tstop, nsize;

	if ( !PadicRuntimeActive ) return(1);
	aux = GetPadicAux();
	if ( aux == 0 ) return(1);

	t = term + *term;
	nsize = ABS(t[-1]);
	tstop = t - nsize;

	if ( nsize == 3 && tstop[0] == 1 && tstop[1] == 1 ) {
		t = term + 1;
		while ( t < tstop ) {
			if ( *t == PADICFUN && t + t[1] == tstop && TestPadic(t) ) {
				return(Generator(BHEAD term,level));
			}
			t += t[1];
		}
	}

	FormRatToMpq(aux->q1,(UWORD *)tstop,(t[-1]));
	padic_set_mpq(aux->p1,aux->q1,PadicContext);
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
	#] Compiler/runtime statements :
	#[ Sorting :
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
		if ( size1 < 0 ) padic_neg(aux->p1,aux->p1,PadicContext);
		UnpackPadic(aux,aux->p2,fun2);
		if ( size2 < 0 ) padic_neg(aux->p2,aux->p2,PadicContext);
	}
	else if ( AT.SortPadicMode == 1 ) {
		fun1 = s1+1; while ( fun1 < coef1 && fun1[0] != PADICFUN ) fun1 += fun1[1];
		UnpackPadic(aux,aux->p1,fun1);
		if ( size1 < 0 ) padic_neg(aux->p1,aux->p1,PadicContext);
		FormRatToMpq(aux->q1,(UWORD *)coef2,size2);
		padic_set_mpq(aux->p2,aux->q1,PadicContext);
	}
	else if ( AT.SortPadicMode == 2 ) {
		fun2 = s2+1; while ( fun2 < coef2 && fun2[0] != PADICFUN ) fun2 += fun2[1];
		UnpackPadic(aux,aux->p2,fun2);
		if ( size2 < 0 ) padic_neg(aux->p2,aux->p2,PadicContext);
		FormRatToMpq(aux->q1,(UWORD *)coef1,size1);
		padic_set_mpq(aux->p1,aux->q1,PadicContext);
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
		if ( size1 < 0 ) padic_neg(aux->p1,aux->p1,PadicContext);
		UnpackPadic(aux,aux->p2,fun2);
		if ( size2 < 0 ) padic_neg(aux->p2,aux->p2,PadicContext);
	}
	else if ( AT.SortPadicMode == 1 ) {
		fun1 = term1+1; while ( fun1 < coef1 && fun1[0] != PADICFUN ) fun1 += fun1[1];
		UnpackPadic(aux,aux->p1,fun1);
		if ( size1 < 0 ) padic_neg(aux->p1,aux->p1,PadicContext);
		FormRatToMpq(aux->q1,(UWORD *)coef2,size2);
		padic_set_mpq(aux->p2,aux->q1,PadicContext);
	}
	else if ( AT.SortPadicMode == 2 ) {
		fun2 = term2+1; while ( fun2 < coef2 && fun2[0] != PADICFUN ) fun2 += fun2[1];
		FormRatToMpq(aux->q1,(UWORD *)coef1,size1);
		padic_set_mpq(aux->p1,aux->q1,PadicContext);
		UnpackPadic(aux,aux->p2,fun2);
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
	#] Sorting :
*/
