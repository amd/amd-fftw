/*
 * Copyright (c) 2003, 2007-14 Matteo Frigo
 * Copyright (c) 2003, 2007-14 Massachusetts Institute of Technology
 * Copyright (C) 2019-2022, Advanced Micro Devices, Inc. All Rights Reserved.
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


#include "dft/dft.h"
#ifdef AMD_FAST_PLANNER
#include "simd-support/simd-common.h"
#endif
#include <stddef.h>

static void destroy(problem *ego_)
{
     problem_dft *ego = (problem_dft *) ego_;
     X(tensor_destroy2)(ego->vecsz, ego->sz);
     X(ifree)(ego_);
}

static void hash(const problem *p_, md5 *m)
{
     const problem_dft *p = (const problem_dft *) p_;
     X(md5puts)(m, "dft");
     X(md5int)(m, p->ri == p->ro);
     X(md5INT)(m, p->ii - p->ri);
     X(md5INT)(m, p->io - p->ro);
     X(md5int)(m, X(ialignment_of)(p->ri));
     X(md5int)(m, X(ialignment_of)(p->ii));
     X(md5int)(m, X(ialignment_of)(p->ro));
     X(md5int)(m, X(ialignment_of)(p->io));

#ifdef AMD_FAST_PLANNER
     X(md5int)(m, p->sz->rnk);
     if (FINITE_RNK(p->sz->rnk)) {
          int x1, x2, x3, x4, x5, x6, x7, x8, x9;
	  for (int i = 0; i < p->sz->rnk; ++i) {
	       X(md5INT)(m, p->sz->dims[i].n);
	       x1 = (p->sz->dims[i].is == p->sz->dims[i].os);
	       x2 = (p->sz->dims[i].is == 2);
	       x3 = (p->sz->dims[i].os == 2);
	       x4 = !((p->sz->dims[i].is * sizeof(R)) % ALIGNMENT);
	       x5 = !((p->sz->dims[i].os * sizeof(R)) % ALIGNMENT);
	       x6 = !((p->sz->dims[i].is * sizeof(R)) % ALIGNMENTA);
	       x7 = !((p->sz->dims[i].os * sizeof(R)) % ALIGNMENTA);
#ifdef AMD_FAST_PLANNING_HASH_V1
	       if (i == 0)
	       {
		       if (p->vecsz->rnk > 0)
			x8 = (p->sz->dims[i].is <= p->vecsz->dims[i].is);
		       else
			x8 = (p->sz->dims[i].is <= 0);
	       }
	       else
	        x8 = 0;
	       x9 = (x8<<7) | (x7<<6) | (x6<<5) | (x5<<4) | (x4<<3) | (x3<<2) | (x2<<1) | x1;
#else //AMD_FAST_PLANNING_HASH_V2
	       x9 = (x7<<6) | (x6<<5) | (x5<<4) | (x4<<3) | (x3<<2) | (x2<<1) | x1;
#endif
	       X(md5INT)(m, x9);
	  }
     }
     int max_ind = X(tensor_max_index)(p->sz);
     X(md5int)(m, p->vecsz->rnk);
     if (FINITE_RNK(p->vecsz->rnk)) {
          int x1=0, x2=0, x3=0, x4, x5, x6, x7, x8, x9, x10=0;
	  for (int i = 0; i < p->vecsz->rnk; ++i) {
	       X(md5INT)(m, p->vecsz->dims[i].n);
	       x1 = (p->vecsz->dims[i].is == p->vecsz->dims[i].os);
	       x2 = (p->vecsz->dims[i].is == 2);
	       x3 = (p->vecsz->dims[i].os == 2);
	       if (x1)
	       {
		       x10 = (X(iabs)(p->vecsz->dims[i].is) < max_ind);
	       }
	       x4 = !((p->vecsz->dims[i].is * sizeof(R)) % ALIGNMENT);
	       x5 = !((p->vecsz->dims[i].os * sizeof(R)) % ALIGNMENT);
	       x6 = !((p->vecsz->dims[i].is * sizeof(R)) % ALIGNMENTA);
	       x7 = !((p->vecsz->dims[i].os * sizeof(R)) % ALIGNMENTA);
	       x9 = (x7<<6) | (x6<<5) | (x5<<4) | (x4<<3) | (x3<<2) | (x2<<1) | x1;
	       X(md5INT)(m, x9);
	  }
	       X(md5int)(m, x10);
     }
#else
     X(tensor_md5)(m, p->sz);
     X(tensor_md5)(m, p->vecsz);
#endif
}

static void print(const problem *ego_, printer *p)
{
     const problem_dft *ego = (const problem_dft *) ego_;
     p->print(p, "(dft %d %d %d %D %D %T %T)", 
	      ego->ri == ego->ro,
	      X(ialignment_of)(ego->ri),
	      X(ialignment_of)(ego->ro),
	      (INT)(ego->ii - ego->ri), 
	      (INT)(ego->io - ego->ro),
	      ego->sz,
	      ego->vecsz);
}

static void zero(const problem *ego_)
{
     const problem_dft *ego = (const problem_dft *) ego_;
     tensor *sz = X(tensor_append)(ego->vecsz, ego->sz);
     X(dft_zerotens)(sz, UNTAINT(ego->ri), UNTAINT(ego->ii));
     X(tensor_destroy)(sz);
}

static const problem_adt padt =
{
     PROBLEM_DFT,
     hash,
     zero,
     print,
     destroy
};

problem *X(mkproblem_dft)(const tensor *sz, const tensor *vecsz,
			  R *ri, R *ii, R *ro, R *io)
{
     problem_dft *ego;

     /* enforce pointer equality if untainted pointers are equal */
     if (UNTAINT(ri) == UNTAINT(ro))
	  ri = ro = JOIN_TAINT(ri, ro);
     if (UNTAINT(ii) == UNTAINT(io))
	  ii = io = JOIN_TAINT(ii, io);

     /* more correctness conditions: */
     A(TAINTOF(ri) == TAINTOF(ii));
     A(TAINTOF(ro) == TAINTOF(io));

     A(X(tensor_kosherp)(sz));
     A(X(tensor_kosherp)(vecsz));

     if (ri == ro || ii == io) {
	  /* If either real or imag pointers are in place, both must be. */
	  if (ri != ro || ii != io || !X(tensor_inplace_locations)(sz, vecsz))
	       return X(mkproblem_unsolvable)();
     }

     ego = (problem_dft *)X(mkproblem)(sizeof(problem_dft), &padt);

     ego->sz = X(tensor_compress)(sz);
     ego->vecsz = X(tensor_compress_contiguous)(vecsz);
     ego->ri = ri;
     ego->ii = ii;
     ego->ro = ro;
     ego->io = io;

     A(FINITE_RNK(ego->sz->rnk));
     return &(ego->super);
}

/* Same as X(mkproblem_dft), but also destroy input tensors. */
#ifdef AMD_FMV_AUTO
__attribute__((target_clones(TARGET_STRINGS)))
#endif
problem *X(mkproblem_dft_d)(tensor *sz, tensor *vecsz,
			    R *ri, R *ii, R *ro, R *io)
{
     problem *p = X(mkproblem_dft)(sz, vecsz, ri, ii, ro, io);
     X(tensor_destroy2)(vecsz, sz);
     return p;
}
