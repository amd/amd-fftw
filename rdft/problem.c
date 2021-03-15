/*
 * Copyright (c) 2003, 2007-14 Matteo Frigo
 * Copyright (c) 2003, 2007-14 Massachusetts Institute of Technology
 * Copyright (C) 2019-2020, Advanced Micro Devices, Inc. All Rights Reserved.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */


#include "rdft/rdft.h"
#ifdef AMD_FAST_PLANNER
#include "simd-support/simd-common.h"
#endif
#include <stddef.h>

#ifdef AMD_FAST_PLANNER
#define MAXRRNK 32 /* FIXME: should malloc() */
typedef struct {
     plan_rdft super;
     INT vl;
     int rnk;
     iodim d[MAXRRNK];
     const char *nam;
} P;
#endif

static void destroy(problem *ego_)
{
     problem_rdft *ego = (problem_rdft *) ego_;
#if !defined(STRUCT_HACK_C99) && !defined(STRUCT_HACK_KR)
     X(ifree0)(ego->kind);
#endif
     X(tensor_destroy2)(ego->vecsz, ego->sz);
     X(ifree)(ego_);
}

static void kind_hash(md5 *m, const rdft_kind *kind, int rnk)
{
     int i;
     for (i = 0; i < rnk; ++i)
	  X(md5int)(m, kind[i]);
}

#ifdef AMD_FAST_PLANNER
static int fill_iodim(P *pln, const problem_rdft *p)
{
     int i;
     const tensor *vecsz = p->vecsz;

     pln->vl = 1;
     pln->rnk = 0;
     for (i = 0; i < vecsz->rnk; ++i) {
	  /* extract contiguous dimensions */
	  if (pln->vl == 1 &&
	      vecsz->dims[i].is == 1 && vecsz->dims[i].os == 1) 
	       pln->vl = vecsz->dims[i].n;
	  else if (pln->rnk == MAXRRNK) 
	       return 0;
	  else 
	       pln->d[pln->rnk++] = vecsz->dims[i];
     }

     return 1;
}
#endif

static void hash(const problem *p_, md5 *m)
{
     const problem_rdft *p = (const problem_rdft *) p_;
     X(md5puts)(m, "rdft");
     X(md5int)(m, p->I == p->O);
     kind_hash(m, p->kind, p->sz->rnk);
     X(md5int)(m, X(ialignment_of)(p->I));
     X(md5int)(m, X(ialignment_of)(p->O));
#ifdef AMD_FAST_PLANNER
     X(md5int)(m, p->sz->rnk);
     if (FINITE_RNK(p->sz->rnk)) {
          int x1, x2, x3, x4, x5, x6, x7, x8, x9, x10;
	  if (p->sz->rnk > 0)
	  {
		  x8 = 1;
		  x10 = X(is_prime)(p->sz->dims[0].n);
	  }
	  else
	  {
		  x8 = 0;
		  x10 = 0;
	  }
	  for (int i = 0; i < p->sz->rnk; ++i) {
	       X(md5INT)(m, p->sz->dims[i].n);
	       x1 = (p->sz->dims[i].is == p->sz->dims[i].os);
	       x2 = (p->sz->dims[i].is == 2);
	       x3 = (p->sz->dims[i].os == 2);
	       x4 = !((p->sz->dims[i].is * sizeof(R)) % ALIGNMENT);
	       x5 = !((p->sz->dims[i].os * sizeof(R)) % ALIGNMENT);
	       x6 = !((p->sz->dims[i].is * sizeof(R)) % ALIGNMENTA);
	       x7 = !((p->sz->dims[i].os * sizeof(R)) % ALIGNMENTA);
	       x8 = x8 & ((p->sz->dims[i].is <= 2) && (p->sz->dims[i].os > 2));
	       x9 = (x8<<7) | (x7<<6) | (x6<<5) | (x5<<4) | (x4<<3) | (x3<<2) | (x2<<1) | x1;
	       X(md5INT)(m, x9);
	  }
	  x10 = (x10<<1) | x8;
	  X(md5int)(m, x10);
     }
     X(md5int)(m, p->vecsz->rnk);
     if (FINITE_RNK(p->vecsz->rnk)) {
          int x1, x2, x3, x4, x5, x6, x7, x8, x9;
	  P pln;
	  for (int i = 0; i < p->vecsz->rnk; ++i) {
	       X(md5INT)(m, p->vecsz->dims[i].n);
	       x1 = (p->vecsz->dims[i].is == p->vecsz->dims[i].os);
	       x2 = (p->vecsz->dims[i].is == 2);
	       x3 = (p->vecsz->dims[i].os == 2);
	       x4 = !((p->vecsz->dims[i].is * sizeof(R)) % ALIGNMENT);
	       x5 = !((p->vecsz->dims[i].os * sizeof(R)) % ALIGNMENT);
	       x6 = !((p->vecsz->dims[i].is * sizeof(R)) % ALIGNMENTA);
	       x7 = !((p->vecsz->dims[i].os * sizeof(R)) % ALIGNMENTA);
	       x9 = (x7<<6) | (x6<<5) | (x5<<4) | (x4<<3) | (x3<<2) | (x2<<1) | x1;
	       X(md5INT)(m, x9);
	  }
	  fill_iodim(&pln, p);
	  x1 = (pln.vl > 2);
	  if (pln.rnk >= 2)
	  {
		  int rnk = pln.rnk;
		  x2 = (pln.d[rnk-2].n == pln.d[rnk-1].n &&
	                pln.d[rnk-2].is == pln.d[rnk-1].os &&
                        pln.d[rnk-2].os == pln.d[rnk-1].is);
		  x3 = (X(iabs)(pln.d[rnk-2].is) <= X(iabs)(pln.d[rnk-1].is) ||
			X(iabs)(pln.d[rnk-2].os) <= X(iabs)(pln.d[rnk-1].os));
		  x4 = (X(compute_tilesz)(pln.vl, 1) > 4);
		  x5 = (X(compute_tilesz)(pln.vl, 2) > 4);

	  }
	  else
	  {
		  x2 = 0;
		  x3 = 0;
		  x4 = 0;
		  x5 = 0;
	  }
	  x9 = (x5<<4) | (x4<<3) | (x3<<2) | (x2<<1) | x1;
	  X(md5int)(m, x9);
     }
#else
     X(tensor_md5)(m, p->sz);
     X(tensor_md5)(m, p->vecsz);
#endif
}

static void recur(const iodim *dims, int rnk, R *I)
{
     if (rnk == RNK_MINFTY)
          return;
     else if (rnk == 0)
          I[0] = K(0.0);
     else if (rnk > 0) {
          INT i, n = dims[0].n, is = dims[0].is;

	  if (rnk == 1) {
	       /* this case is redundant but faster */
	       for (i = 0; i < n; ++i)
		    I[i * is] = K(0.0);
	  } else {
	       for (i = 0; i < n; ++i)
		    recur(dims + 1, rnk - 1, I + i * is);
	  }
     }
}

void X(rdft_zerotens)(tensor *sz, R *I)
{
     recur(sz->dims, sz->rnk, I);
}

#define KSTR_LEN 8

const char *X(rdft_kind_str)(rdft_kind kind)
{
     static const char kstr[][KSTR_LEN] = {
	  "r2hc", "r2hc01", "r2hc10", "r2hc11",
	  "hc2r", "hc2r01", "hc2r10", "hc2r11",
	  "dht",
	  "redft00", "redft01", "redft10", "redft11",
	  "rodft00", "rodft01", "rodft10", "rodft11"
     };
     A(kind >= 0 && kind < sizeof(kstr) / KSTR_LEN);
     return kstr[kind];
}

static void print(const problem *ego_, printer *p)
{
     const problem_rdft *ego = (const problem_rdft *) ego_;
     int i;
     p->print(p, "(rdft %d %D %T %T", 
	      X(ialignment_of)(ego->I),
	      (INT)(ego->O - ego->I), 
	      ego->sz,
	      ego->vecsz);
     for (i = 0; i < ego->sz->rnk; ++i)
	  p->print(p, " %d", (int)ego->kind[i]);
     p->print(p, ")");
}

static void zero(const problem *ego_)
{
     const problem_rdft *ego = (const problem_rdft *) ego_;
     tensor *sz = X(tensor_append)(ego->vecsz, ego->sz);
     X(rdft_zerotens)(sz, UNTAINT(ego->I));
     X(tensor_destroy)(sz);
}

static const problem_adt padt =
{
     PROBLEM_RDFT,
     hash,
     zero,
     print,
     destroy
};

/* Dimensions of size 1 that are not REDFT/RODFT are no-ops and can be
   eliminated.  REDFT/RODFT unit dimensions often have factors of 2.0
   and suchlike from normalization and phases, although in principle
   these constant factors from different dimensions could be combined. */
static int nontrivial(const iodim *d, rdft_kind kind)
{
     return (d->n > 1 || kind == R2HC11 || kind == HC2R11
	     || (REODFT_KINDP(kind) && kind != REDFT01 && kind != RODFT01));
}

problem *X(mkproblem_rdft)(const tensor *sz, const tensor *vecsz,
			   R *I, R *O, const rdft_kind *kind)
{
     problem_rdft *ego;
     int rnk = sz->rnk;
     int i;

     A(X(tensor_kosherp)(sz));
     A(X(tensor_kosherp)(vecsz));
     A(FINITE_RNK(sz->rnk));

     if (UNTAINT(I) == UNTAINT(O))
	  I = O = JOIN_TAINT(I, O);

     if (I == O && !X(tensor_inplace_locations)(sz, vecsz))
	  return X(mkproblem_unsolvable)();

     for (i = rnk = 0; i < sz->rnk; ++i) {
          A(sz->dims[i].n > 0);
          if (nontrivial(sz->dims + i, kind[i]))
               ++rnk;
     }

#if defined(STRUCT_HACK_KR)
     ego = (problem_rdft *) X(mkproblem)(sizeof(problem_rdft)
					 + sizeof(rdft_kind)
					 * (rnk > 0 ? rnk - 1u : 0u), &padt);
#elif defined(STRUCT_HACK_C99)
     ego = (problem_rdft *) X(mkproblem)(sizeof(problem_rdft)
					 + sizeof(rdft_kind) * (unsigned)rnk, &padt);
#else
     ego = (problem_rdft *) X(mkproblem)(sizeof(problem_rdft), &padt);
     ego->kind = (rdft_kind *) MALLOC(sizeof(rdft_kind) * (unsigned)rnk, PROBLEMS);
#endif

     /* do compression and sorting as in X(tensor_compress), but take
	transform kind into account (sigh) */
     ego->sz = X(mktensor)(rnk);
     for (i = rnk = 0; i < sz->rnk; ++i) {
          if (nontrivial(sz->dims + i, kind[i])) {
	       ego->kind[rnk] = kind[i];
               ego->sz->dims[rnk++] = sz->dims[i];
	  }
     }
     for (i = 0; i + 1 < rnk; ++i) {
	  int j;
	  for (j = i + 1; j < rnk; ++j)
	       if (X(dimcmp)(ego->sz->dims + i, ego->sz->dims + j) > 0) {
		    iodim dswap;
		    rdft_kind kswap;
		    dswap = ego->sz->dims[i];
		    ego->sz->dims[i] = ego->sz->dims[j];
		    ego->sz->dims[j] = dswap;
		    kswap = ego->kind[i];
		    ego->kind[i] = ego->kind[j];
		    ego->kind[j] = kswap;
	       }
     }

     for (i = 0; i < rnk; ++i)
	  if (ego->sz->dims[i].n == 2 && (ego->kind[i] == REDFT00
					  || ego->kind[i] == DHT
					  || ego->kind[i] == HC2R))
	       ego->kind[i] = R2HC; /* size-2 transforms are equivalent */

     ego->vecsz = X(tensor_compress_contiguous)(vecsz);
     ego->I = I;
     ego->O = O;

     A(FINITE_RNK(ego->sz->rnk));

     return &(ego->super);
}

/* Same as X(mkproblem_rdft), but also destroy input tensors. */
problem *X(mkproblem_rdft_d)(tensor *sz, tensor *vecsz,
			     R *I, R *O, const rdft_kind *kind)
{
     problem *p = X(mkproblem_rdft)(sz, vecsz, I, O, kind);
     X(tensor_destroy2)(vecsz, sz);
     return p;
}

/* As above, but for rnk <= 1 only and takes a scalar kind parameter */
problem *X(mkproblem_rdft_1)(const tensor *sz, const tensor *vecsz,
			     R *I, R *O, rdft_kind kind)
{
     A(sz->rnk <= 1);
     return X(mkproblem_rdft)(sz, vecsz, I, O, &kind);
}

problem *X(mkproblem_rdft_1_d)(tensor *sz, tensor *vecsz,
			       R *I, R *O, rdft_kind kind)
{
     A(sz->rnk <= 1);
     return X(mkproblem_rdft_d)(sz, vecsz, I, O, &kind);
}

/* create a zero-dimensional problem */
problem *X(mkproblem_rdft_0_d)(tensor *vecsz, R *I, R *O)
{
     return X(mkproblem_rdft_d)(X(mktensor_0d)(), vecsz, I, O, 
				(const rdft_kind *)0);
}
