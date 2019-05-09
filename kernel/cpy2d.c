/*
 * Copyright (c) 2003, 2007-14 Matteo Frigo
 * Copyright (c) 2003, 2007-14 Massachusetts Institute of Technology
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

/* out of place 2D copy routines */
#include "kernel/ifftw.h"

#ifdef AMD_OPT_ALL
#include "immintrin.h"
#endif

#if defined(__x86_64__) || defined(_M_X64) || defined(_M_AMD64)
#  ifdef HAVE_XMMINTRIN_H
#    include <xmmintrin.h>
#    define WIDE_TYPE __m128
#  endif
#endif

#ifndef WIDE_TYPE
/* fall back to double, which means that WIDE_TYPE will be unused */
#  define WIDE_TYPE double
#endif

#ifdef AMD_OPT_ALL //AMD optimized routines

#ifdef AMD_OPT_IN_PLACE_1D_CPY2D_STABLE_C
void X(cpy2d)(R *I, R *O,
	      INT n0, INT is0, INT os0,
	      INT n1, INT is1, INT os1,
	      INT vl)
{
     INT i0, i1, v;

     switch (vl) {
	 case 1:
	      for (i1 = 0; i1 < n1; ++i1)
		   for (i0 = 0; i0 < n0; ++i0) {
			R x0 = I[i0 * is0 + i1 * is1];
			O[i0 * os0 + i1 * os1] = x0;
		   }
	      break;
	 case 2:
	      if (1
		  && (2 * sizeof(R) == sizeof(WIDE_TYPE))
		  && (sizeof(WIDE_TYPE) > sizeof(double))
		  && (((size_t)I) % sizeof(WIDE_TYPE) == 0)
		  && (((size_t)O) % sizeof(WIDE_TYPE) == 0)
		  && ((is0 & 1) == 0)
		  && ((is1 & 1) == 0)
		  && ((os0 & 1) == 0)
		  && ((os1 & 1) == 0)) {
		   /* copy R[2] as WIDE_TYPE if WIDE_TYPE is large
		      enough to hold R[2], and if the input is
		      properly aligned.  This is a win when R==double
		      and WIDE_TYPE is 128 bits. */
		   for (i1 = 0; i1 < n1; ++i1)
			for (i0 = 0; i0 < n0; ++i0) {
			     *(WIDE_TYPE *)&O[i0 * os0 + i1 * os1] =
				  *(WIDE_TYPE *)&I[i0 * is0 + i1 * is1];
			}
	      } else if (1
		  && (2 * sizeof(R) == sizeof(double))
		  && (((size_t)I) % sizeof(double) == 0)
		  && (((size_t)O) % sizeof(double) == 0)
		  && ((is0 & 1) == 0)
		  && ((is1 & 1) == 0)
		  && ((os0 & 1) == 0)
		  && ((os1 & 1) == 0)) {
		   /* copy R[2] as double if double is large enough to
		      hold R[2], and if the input is properly aligned.
		      This case applies when R==float */
		   for (i1 = 0; i1 < n1; ++i1)
			for (i0 = 0; i0 < n0; ++i0) {
			     *(double *)&O[i0 * os0 + i1 * os1] =
				  *(double *)&I[i0 * is0 + i1 * is1];
			}
	      } else {
		      INT t0, t1, t2, t3;
		      INT n0_rem = n0&0x1, n1_rem = n1&0x1;
		      t0 = (is0==2);
		      t1 = (os0==2);
		      t2 = (is1==2);
		      t3 = (os1==2);
		      
		      switch(t0 | (t1 << 1) | (t2 << 2) | (t3 << 3))
		      {
			      case 1://only is0 is 2. 256-bit contiguous read possible
			      n0 = n0 - n0_rem;
			      for (i1 = 0; i1 < n1; ++i1) {
				      for (i0 = 0; i0 < n0; i0+=2) {
					      R x0 = I[i0 * is0 + i1 * is1];
					      R x1 = I[i0 * is0 + i1 * is1 + 1];
					      R x2 = I[i0 * is0 + i1 * is1 + 2];
					      R x3 = I[i0 * is0 + i1 * is1 + 3];

					      O[i0 * os0 + i1 * os1] = x0;
					      O[i0 * os0 + i1 * os1 + 1] = x1;
					      O[(i0+1) * os0 + i1 * os1] = x2;
					      O[(i0+1) * os0 + i1 * os1 + 1] = x3;
				      }
				      if (n0_rem)
				      {
					      R x0 = I[i0 * is0 + i1 * is1];
					      R x1 = I[i0 * is0 + i1 * is1 + 1];
					      O[i0 * os0 + i1 * os1] = x0;
					      O[i0 * os0 + i1 * os1 + 1] = x1;
				      }
			      }
			      break;

			      case 2://only os0 is 2. 256-bit contiguous write possible
			      n0 = n0 - n0_rem;
			      for (i1 = 0; i1 < n1; ++i1) {
				      for (i0 = 0; i0 < n0; i0+=2) {
					      R x0 = I[i0 * is0 + i1 * is1];
					      R x1 = I[i0 * is0 + i1 * is1 + 1];
					      R x2 = I[(i0+1) * is0 + i1 * is1];
					      R x3 = I[(i0+1) * is0 + i1 * is1 + 1];

					      O[i0 * os0 + i1 * os1] = x0;
					      O[i0 * os0 + i1 * os1 + 1] = x1;
					      O[i0 * os0 + i1 * os1 + 2] = x2;
					      O[i0 * os0 + i1 * os1 + 3] = x3;
				      }
				      if (n0_rem)
				      {
					      R x0 = I[i0 * is0 + i1 * is1];
					      R x1 = I[i0 * is0 + i1 * is1 + 1];
					      O[i0 * os0 + i1 * os1] = x0;
					      O[i0 * os0 + i1 * os1 + 1] = x1;
				      }
			      }
			      break;

			      case 6://os0=2 and is1=2. Both 256-bit read and 256-bit write possible
			      n0 = n0 - n0_rem;
			      n1 = n1 - n1_rem;
			      for (i1 = 0; i1 < n1; i1+=2)
			      {
				      for (i0 = 0; i0 < n0; i0+=2) {
					      R x0 = I[i0 * is0 + i1 * is1];
					      R x1 = I[i0 * is0 + i1 * is1 + 1];
					      R x2 = I[i0 * is0 + (i1+1) * is1];
					      R x3 = I[i0 * is0 + (i1+1) * is1 + 1];
					      R x4 = I[(i0+1) * is0 + i1 * is1];
					      R x5 = I[(i0+1) * is0 + i1 * is1 + 1];
					      R x6 = I[(i0+1) * is0 + (i1+1) * is1];
					      R x7 = I[(i0+1) * is0 + (i1+1) * is1 + 1];
					      O[i0 * os0 + i1 * os1] = x0;
					      O[i0 * os0 + i1 * os1 + 1] = x1;
					      O[(i0+1) * os0 + i1 * os1] = x4;
					      O[(i0+1) * os0 + i1 * os1 + 1] = x5;
					      O[i0 * os0 + (i1+1) * os1] = x2;
					      O[i0 * os0 + (i1+1) * os1 + 1] = x3;
					      O[(i0+1) * os0 + (i1+1) * os1] = x6;
					      O[(i0+1) * os0 + (i1+1) * os1 + 1] = x7;
				      }
				      if (n0_rem)
				      {
					      R x0 = I[i0 * is0 + i1 * is1];
					      R x1 = I[i0 * is0 + i1 * is1 + 1];
					      R x2 = I[i0 * is0 + (i1+1) * is1];
					      R x3 = I[i0 * is0 + (i1+1) * is1 + 1];
					      O[i0 * os0 + i1 * os1] = x0;
					      O[i0 * os0 + i1 * os1 + 1] = x1;
					      O[i0 * os0 + (i1+1) * os1] = x2;
					      O[i0 * os0 + (i1+1) * os1 + 1] = x3;
				      }
			      }
			      n0 += n0_rem;
			      if (n1_rem)
			      {
				      for (i0 = 0; i0 < n0; ++i0) {
					      R x0 = I[i0 * is0 + i1 * is1];
					      R x1 = I[i0 * is0 + i1 * is1 + 1];
					      O[i0 * os0 + i1 * os1] = x0;
					      O[i0 * os0 + i1 * os1 + 1] = x1;
				      }
			      }
			      break;

			      case 8://only os1 is 2. 256-bit contiguous write possible
			      n1 = n1 - n1_rem;
			      for (i1 = 0; i1 < n1; i1+=2) {
				      for (i0 = 0; i0 < n0; ++i0) {
					      R x0 = I[i0 * is0 + i1 * is1];
					      R x1 = I[i0 * is0 + i1 * is1 + 1];
					      R x2 = I[i0 * is0 + (i1+1) * is1];
					      R x3 = I[i0 * is0 + (i1+1) * is1 + 1];

					      O[i0 * os0 + i1 * os1] = x0;
					      O[i0 * os0 + i1 * os1 + 1] = x1;
					      O[i0 * os0 + i1 * os1 + 2] = x2;
					      O[i0 * os0 + i1 * os1 + 3] = x3;
				      }
			      }
			      if (n1_rem)
			      {
				      for (i0 = 0; i0 < n0; ++i0) {
					      R x0 = I[i0 * is0 + i1 * is1];
					      R x1 = I[i0 * is0 + i1 * is1 + 1];
					      O[i0 * os0 + i1 * os1] = x0;
					      O[i0 * os0 + i1 * os1 + 1] = x1;
				      }
			      }
			      break;
			      
			      case 9://is0=2 and os1=2. Both 256-bit read and 256-bit write possible
			      n0 = n0 - n0_rem;
			      n1 = n1 - n1_rem;
			      for (i1 = 0; i1 < n1; i1+=2)
			      {
				      for (i0 = 0; i0 < n0; i0+=2) {
					      R x0 = I[i0 * is0 + i1 * is1];
					      R x1 = I[i0 * is0 + i1 * is1 + 1];
					      R x2 = I[(i0+1) * is0 + i1 * is1];
					      R x3 = I[(i0+1) * is0 + i1 * is1 + 1];
					      R x4 = I[i0 * is0 + (i1+1) * is1];
					      R x5 = I[i0 * is0 + (i1+1) * is1 + 1];
					      R x6 = I[(i0+1) * is0 + (i1+1) * is1];
					      R x7 = I[(i0+1) * is0 + (i1+1) * is1 + 1];
					      O[i0 * os0 + i1 * os1] = x0;
					      O[i0 * os0 + i1 * os1 + 1] = x1;
					      O[i0 * os0 + (i1+1) * os1] = x4;
					      O[i0 * os0 + (i1+1) * os1 + 1] = x5;
					      O[(i0+1) * os0 + i1 * os1] = x2;
					      O[(i0+1) * os0 + i1 * os1 + 1] = x3;
					      O[(i0+1) * os0 + (i1+1) * os1] = x6;
					      O[(i0+1) * os0 + (i1+1) * os1 + 1] = x7;
				      }
				      if (n0_rem)
				      {
					      R x0 = I[i0 * is0 + i1 * is1];
					      R x1 = I[i0 * is0 + i1 * is1 + 1];
					      R x2 = I[i0 * is0 + (i1+1) * is1];
					      R x3 = I[i0 * is0 + (i1+1) * is1 + 1];
					      O[i0 * os0 + i1 * os1] = x0;
					      O[i0 * os0 + i1 * os1 + 1] = x1;
					      O[i0 * os0 + (i1+1) * os1] = x2;
					      O[i0 * os0 + (i1+1) * os1 + 1] = x3;
				      }
			      }
			      if (n1_rem)
			      {
				      n0 += n0_rem;
				      for (i0 = 0; i0 < n0; ++i0) {
					      R x0 = I[i0 * is0 + i1 * is1];
					      R x1 = I[i0 * is0 + i1 * is1 + 1];
					      O[i0 * os0 + i1 * os1] = x0;
					      O[i0 * os0 + i1 * os1 + 1] = x1;
				      }
			      }
			      break;

			      case 3://is0=2 and os0=2. Both 256-bit read and 256-bit write possible
			      case 7://is0=2 and os0=2. Also is1=2. Both 256-bit read and 256-bit write possible
			      case 11://is0=2 and os0=2. Also os1=2. Both 256-bit read and 256-bit write possible
			      case 15://is0=2 and os0=2. Also is1=2, os1=2. Both 256-bit read and 256-bit write possible
			      n0 = n0 - n0_rem;
			      for (i1 = 0; i1 < n1; ++i1)
			      {
				      for (i0 = 0; i0 < n0; i0+=2) {
					      R x0 = I[i0 * is0 + i1 * is1];
					      R x1 = I[i0 * is0 + i1 * is1 + 1];
					      R x2 = I[(i0+1) * is0 + i1 * is1];
					      R x3 = I[(i0+1) * is0 + i1 * is1 + 1];
					      O[i0 * os0 + i1 * os1] = x0;
					      O[i0 * os0 + i1 * os1 + 1] = x1;
					      O[(i0+1) * os0 + i1 * os1] = x2;
					      O[(i0+1) * os0 + i1 * os1 + 1] = x3;
				      }
				      if (n0_rem)
				      {
					      R x0 = I[i0 * is0 + i1 * is1];
					      R x1 = I[i0 * is0 + i1 * is1 + 1];
					      O[i0 * os0 + i1 * os1] = x0;
					      O[i0 * os0 + i1 * os1 + 1] = x1;
				      }
			      }
			      break;
			      
			      default:
			      for (i1 = 0; i1 < n1; ++i1)
				      for (i0 = 0; i0 < n0; ++i0) {
					      R x0 = I[i0 * is0 + i1 * is1];
					      R x1 = I[i0 * is0 + i1 * is1 + 1];
					      O[i0 * os0 + i1 * os1] = x0;
					      O[i0 * os0 + i1 * os1 + 1] = x1;
				      }
			      break;
		      }
	      }
	      break;
	 default:
	      for (i1 = 0; i1 < n1; ++i1)
		   for (i0 = 0; i0 < n0; ++i0)
			for (v = 0; v < vl; ++v) {
			     R x0 = I[i0 * is0 + i1 * is1 + v];
			     O[i0 * os0 + i1 * os1 + v] = x0;
			}
	      break;
     }
}
#elif def AMD_OPT_IN_PLACE_1D_CPY2D_STABLE_INTRIN
void X(cpy2d)(R *I, R *O,
	      INT n0, INT is0, INT os0,
	      INT n1, INT is1, INT os1,
	      INT vl)
{
     INT i0, i1, v;

     switch (vl) {
	 case 1:
	      for (i1 = 0; i1 < n1; ++i1)
		   for (i0 = 0; i0 < n0; ++i0) {
			R x0 = I[i0 * is0 + i1 * is1];
			O[i0 * os0 + i1 * os1] = x0;
		   }
	      break;
	 case 2:
	      if (1
		  && (2 * sizeof(R) == sizeof(WIDE_TYPE))
		  && (sizeof(WIDE_TYPE) > sizeof(double))
		  && (((size_t)I) % sizeof(WIDE_TYPE) == 0)
		  && (((size_t)O) % sizeof(WIDE_TYPE) == 0)
		  && ((is0 & 1) == 0)
		  && ((is1 & 1) == 0)
		  && ((os0 & 1) == 0)
		  && ((os1 & 1) == 0)) {
		   /* copy R[2] as WIDE_TYPE if WIDE_TYPE is large
		      enough to hold R[2], and if the input is
		      properly aligned.  This is a win when R==double
		      and WIDE_TYPE is 128 bits. */
		   for (i1 = 0; i1 < n1; ++i1)
			for (i0 = 0; i0 < n0; ++i0) {
			     *(WIDE_TYPE *)&O[i0 * os0 + i1 * os1] =
				  *(WIDE_TYPE *)&I[i0 * is0 + i1 * is1];
			}
	      } else if (1
		  && (2 * sizeof(R) == sizeof(double))
		  && (((size_t)I) % sizeof(double) == 0)
		  && (((size_t)O) % sizeof(double) == 0)
		  && ((is0 & 1) == 0)
		  && ((is1 & 1) == 0)
		  && ((os0 & 1) == 0)
		  && ((os1 & 1) == 0)) {
		   /* copy R[2] as double if double is large enough to
		      hold R[2], and if the input is properly aligned.
		      This case applies when R==float */
		   for (i1 = 0; i1 < n1; ++i1)
			for (i0 = 0; i0 < n0; ++i0) {
			     *(double *)&O[i0 * os0 + i1 * os1] =
				  *(double *)&I[i0 * is0 + i1 * is1];
			}
	      } else {
		      __m256d in1, in2, in3, in4;
		      __m256d out1, out2;
		      __m128d in1_128, in2_128;
		      INT t0, t1, t2, t3;
		      INT n0_rem = n0&0x1, n1_rem = n1&0x1;
		      t0 = (is0==2);
		      t1 = (os0==2);
		      t2 = (is1==2);
		      t3 = (os1==2);
		      
		      switch(t0 | (t1 << 1) | (t2 << 2) | (t3 << 3))
		      {
			      case 1://only is0 is 2. 256-bit contiguous read possible
			      n0 = n0 - n0_rem;
			      for (i1 = 0; i1 < n1; ++i1) {
				      for (i0 = 0; i0 < n0; i0+=2) {
					      in1 = _mm256_loadu_pd((double const *)&I[i0 * is0 + i1 * is1]);
					      in1_128 = _mm256_castpd256_pd128(in1);
					      in2_128 = _mm256_extractf128_pd(in1, 0x1);
					      _mm_storeu_pd((double *)&O[i0 * os0 + i1 * os1], in1_128);
					      _mm_storeu_pd((double *)&O[(i0+1) * os0 + i1 * os1], in2_128);
				      }
				      if (n0_rem)
				      {
					      R x0 = I[i0 * is0 + i1 * is1];
					      R x1 = I[i0 * is0 + i1 * is1 + 1];
					      O[i0 * os0 + i1 * os1] = x0;
					      O[i0 * os0 + i1 * os1 + 1] = x1;
				      }
			      }
			      break;

			      case 2://only os0 is 2. 256-bit contiguous write possible
			      n0 = n0 - n0_rem;
			      for (i1 = 0; i1 < n1; ++i1) {
				      for (i0 = 0; i0 < n0; i0+=2) {
					      in1_128 = _mm_loadu_pd((double const *)&I[i0 * is0 + i1 * is1]);
					      in2_128 = _mm_loadu_pd((double const *)&I[(i0+1) * is0 + i1 * is1]);
					      in1 = _mm256_castpd128_pd256(in1_128);
					      //in2 = _mm256_castpd128_pd256(in2_128);
					      //out1 = _mm256_permute2f128_pd(in1, in2, 0x20);
					      out1 = _mm256_insertf128_pd(in1, in2_128, 1);
					      _mm256_storeu_pd((double *)&O[i0 * os0 + i1 * os1], out1);
				      }
				      if (n0_rem)
				      {
					      R x0 = I[i0 * is0 + i1 * is1];
					      R x1 = I[i0 * is0 + i1 * is1 + 1];
					      O[i0 * os0 + i1 * os1] = x0;
					      O[i0 * os0 + i1 * os1 + 1] = x1;
				      }
			      }
			      break;

			      case 6://os0=2 and is1=2. Both 256-bit read and 256-bit write possible
			      n0 = n0 - n0_rem;
			      n1 = n1 - n1_rem;
			      for (i1 = 0; i1 < n1; i1+=2)
			      {
				      for (i0 = 0; i0 < n0; i0+=2) {
					      in1 = _mm256_loadu_pd((double const *)&I[i0 * is0 + i1 * is1]);
					      in2 = _mm256_loadu_pd((double const *)&I[(i0+1) * is0 + i1 * is1]);

					      //out1 = _mm256_shuffle_pd(in1, in2, 0x33);
					      //out2 = _mm256_shuffle_pd(in1, in2, 0x11);
					      out1 = _mm256_permute2f128_pd(in1, in2, 0x20);
					      out2 = _mm256_permute2f128_pd(in1, in2, 0x31);
					      _mm256_storeu_pd((double *)&O[i0 * os0 + i1 * os1], out1);
					      _mm256_storeu_pd((double *)&O[i0 * os0 + (i1+1) * os1], out2);
				      }
				      if (n0_rem)
				      {
					      R x0 = I[i0 * is0 + i1 * is1];
					      R x1 = I[i0 * is0 + i1 * is1 + 1];
					      R x2 = I[i0 * is0 + (i1+1) * is1];
					      R x3 = I[i0 * is0 + (i1+1) * is1 + 1];
					      O[i0 * os0 + i1 * os1] = x0;
					      O[i0 * os0 + i1 * os1 + 1] = x1;
					      O[i0 * os0 + (i1+1) * os1] = x2;
					      O[i0 * os0 + (i1+1) * os1 + 1] = x3;
				      }
			      }
			      n0 += n0_rem;
			      if (n1_rem)
			      {
				      for (i0 = 0; i0 < n0; ++i0) {
					      R x0 = I[i0 * is0 + i1 * is1];
					      R x1 = I[i0 * is0 + i1 * is1 + 1];
					      O[i0 * os0 + i1 * os1] = x0;
					      O[i0 * os0 + i1 * os1 + 1] = x1;
				      }
			      }
			      break;

			      case 8://only os1 is 2. 256-bit contiguous write possible
			      n1 = n1 - n1_rem;
			      for (i1 = 0; i1 < n1; i1+=2) {
				      for (i0 = 0; i0 < n0; ++i0) {
					      in1_128 = _mm_loadu_pd((double const *)&I[i0 * is0 + i1 * is1]);
					      in2_128 = _mm_loadu_pd((double const *)&I[i0 * is0 + (i1+1) * is1]);
					      in1 = _mm256_castpd128_pd256(in1_128);
					      //in2 = _mm256_castpd128_pd256(in2_128);
					      //out1 = _mm256_permute2f128_pd(in1, in2, 0x20);
					      out1 = _mm256_insertf128_pd(in1, in2_128, 1);
					      _mm256_storeu_pd((double *)&O[i0 * os0 + i1 * os1], out1);
				      }
			      }
			      if (n1_rem)
			      {
				      for (i0 = 0; i0 < n0; ++i0) {
					      R x0 = I[i0 * is0 + i1 * is1];
					      R x1 = I[i0 * is0 + i1 * is1 + 1];
					      O[i0 * os0 + i1 * os1] = x0;
					      O[i0 * os0 + i1 * os1 + 1] = x1;
				      }
			      }
			      break;
			      
			      case 9://is0=2 and os1=2. Both 256-bit read and 256-bit write possible
			      n0 = n0 - n0_rem;
			      n1 = n1 - n1_rem;
			      for (i1 = 0; i1 < n1; i1+=2)
			      {
				      for (i0 = 0; i0 < n0; i0+=2) {
					      in1 = _mm256_loadu_pd((double const *)&I[i0 * is0 + i1 * is1]);
					      in2 = _mm256_loadu_pd((double const *)&I[i0 * is0 + (i1+1) * is1]);

					      //out1 = _mm256_shuffle_pd(in1, in2, 0x33);
					      //out2 = _mm256_shuffle_pd(in1, in2, 0x11);
					      out1 = _mm256_permute2f128_pd(in1, in2, 0x20);
					      out2 = _mm256_permute2f128_pd(in1, in2, 0x31);
					      _mm256_storeu_pd((double *)&O[i0 * os0 + i1 * os1], out1);
					      _mm256_storeu_pd((double *)&O[(i0+1) * os0 + i1 * os1], out2);
				      }
				      if (n0_rem)
				      {
					      R x0 = I[i0 * is0 + i1 * is1];
					      R x1 = I[i0 * is0 + i1 * is1 + 1];
					      R x2 = I[i0 * is0 + (i1+1) * is1];
					      R x3 = I[i0 * is0 + (i1+1) * is1 + 1];
					      O[i0 * os0 + i1 * os1] = x0;
					      O[i0 * os0 + i1 * os1 + 1] = x1;
					      O[i0 * os0 + (i1+1) * os1] = x2;
					      O[i0 * os0 + (i1+1) * os1 + 1] = x3;
				      }
			      }
			      if (n1_rem)
			      {
				      n0 += n0_rem;
				      for (i0 = 0; i0 < n0; ++i0) {
					      R x0 = I[i0 * is0 + i1 * is1];
					      R x1 = I[i0 * is0 + i1 * is1 + 1];
					      O[i0 * os0 + i1 * os1] = x0;
					      O[i0 * os0 + i1 * os1 + 1] = x1;
				      }
			      }
			      break;

			      case 3://is0=2 and os0=2. Both 256-bit read and 256-bit write possible
			      case 7://is0=2 and os0=2. Also is1=2. Both 256-bit read and 256-bit write possible
			      case 11://is0=2 and os0=2. Also os1=2. Both 256-bit read and 256-bit write possible
			      case 15://is0=2 and os0=2. Also is1=2, os1=2. Both 256-bit read and 256-bit write possible
			      n0 = n0 - n0_rem;
			      for (i1 = 0; i1 < n1; ++i1)
			      {
				      for (i0 = 0; i0 < n0; i0+=2) {
					      in1 = _mm256_loadu_pd((double const *)&I[i0 * is0 + i1 * is1]);
					      _mm256_storeu_pd((double *)&O[i0 * os0 + i1 * os1], in1);
				      }
				      if (n0_rem)
				      {
					      R x0 = I[i0 * is0 + i1 * is1];
					      R x1 = I[i0 * is0 + i1 * is1 + 1];
					      O[i0 * os0 + i1 * os1] = x0;
					      O[i0 * os0 + i1 * os1 + 1] = x1;
				      }
			      }
			      break;
			      
			      default:
			      for (i1 = 0; i1 < n1; ++i1)
				      for (i0 = 0; i0 < n0; ++i0) {
					      R x0 = I[i0 * is0 + i1 * is1];
					      R x1 = I[i0 * is0 + i1 * is1 + 1];
					      O[i0 * os0 + i1 * os1] = x0;
					      O[i0 * os0 + i1 * os1 + 1] = x1;
				      }
			      break;
		      }
	      }
	      break;
	 default:
	      for (i1 = 0; i1 < n1; ++i1)
		   for (i0 = 0; i0 < n0; ++i0)
			for (v = 0; v < vl; ++v) {
			     R x0 = I[i0 * is0 + i1 * is1 + v];
			     O[i0 * os0 + i1 * os1 + v] = x0;
			}
	      break;
     }
}
#else //AMD_OPT_IN_PLACE_1D_CPY2D_EXPERIMENTAL_ASM
void X(cpy2d)(R *I, R *O,
	      INT n0, INT is0, INT os0,
	      INT n1, INT is1, INT os1,
	      INT vl)
{
     INT i0, i1, v;

     switch (vl) {
	 case 1:
	      for (i1 = 0; i1 < n1; ++i1)
		   for (i0 = 0; i0 < n0; ++i0) {
			R x0 = I[i0 * is0 + i1 * is1];
			O[i0 * os0 + i1 * os1] = x0;
		   }
	      break;
	 case 2:
	      if (1
		  && (2 * sizeof(R) == sizeof(WIDE_TYPE))
		  && (sizeof(WIDE_TYPE) > sizeof(double))
		  && (((size_t)I) % sizeof(WIDE_TYPE) == 0)
		  && (((size_t)O) % sizeof(WIDE_TYPE) == 0)
		  && ((is0 & 1) == 0)
		  && ((is1 & 1) == 0)
		  && ((os0 & 1) == 0)
		  && ((os1 & 1) == 0)) {
		   /* copy R[2] as WIDE_TYPE if WIDE_TYPE is large
		      enough to hold R[2], and if the input is
		      properly aligned.  This is a win when R==double
		      and WIDE_TYPE is 128 bits. */
		   for (i1 = 0; i1 < n1; ++i1)
			for (i0 = 0; i0 < n0; ++i0) {
			     *(WIDE_TYPE *)&O[i0 * os0 + i1 * os1] =
				  *(WIDE_TYPE *)&I[i0 * is0 + i1 * is1];
			}
	      } else if (1
		  && (2 * sizeof(R) == sizeof(double))
		  && (((size_t)I) % sizeof(double) == 0)
		  && (((size_t)O) % sizeof(double) == 0)
		  && ((is0 & 1) == 0)
		  && ((is1 & 1) == 0)
		  && ((os0 & 1) == 0)
		  && ((os1 & 1) == 0)) {
		   /* copy R[2] as double if double is large enough to
		      hold R[2], and if the input is properly aligned.
		      This case applies when R==float */
		   for (i1 = 0; i1 < n1; ++i1)
			for (i0 = 0; i0 < n0; ++i0) {
			     *(double *)&O[i0 * os0 + i1 * os1] =
				  *(double *)&I[i0 * is0 + i1 * is1];
			}
	      } else {
		      INT t0, t1, t2, t3;
		      INT n0_rem = n0&0x1, n1_rem = n1&0x1;
		      __m256d in1, in2, in3, in4;
		      __m256d out1, out2;
		      __m128d in1_128, in2_128;
		      t0 = (is0==2);
		      t1 = (os0==2);
		      t2 = (is1==2);
		      t3 = (os1==2);
		      
		      switch(t0 | (t1 << 1) | (t2 << 2) | (t3 << 3))
		      {
			      case 1://only is0 is 2. 256-bit contiguous read possible
			      n0 = n0 - n0_rem;
			      _mm256_zeroupper();
			      for (i1 = 0; i1 < n1; ++i1) {
				      t1 = i1 * is1;
				      t2 = i1 * os1;
				      for (i0 = 0; i0 < n0; i0+=2) {
					      __asm__ volatile (
							      //"vzeroall \n\t"
							      "mov %0, %%rcx \n\t"
							      "mov %1, %%rax \n\t"
							      "mul %%rcx \n\t"
							      "mov %2, %%rdx \n\t"
							      "mov %5, %%rsi \n\t"
							      "add %3, %%rsi\n\t"
							      "mov %6, %%rdi\n\t"
							      "add %4, %%rdi\n\t"
							      "vmovupd (%%rsi, %%rax), %%ymm0\n\t"
							      "mov %%rcx, %%rax \n\t"
							      "mul %%rdx \n\t"
							      "vextractf128 $1, %%ymm0, %%xmm1\n\t"
							      "vmovupd %%ymm0, (%%rdi, %%rax)\n\t"
							      "add %%rdx, %%rax\n\t"
							      "vmovupd %%ymm1, (%%rdi, %%rax)\n\t"
							      :// output operands (none)
							      :// input operands
							      "m" (i0), //0
							      "m" (is0),//1
							      "m" (os0),//2
							      "m" (t1), //3
							      "m" (t2), //4
							      "m" (I),  //5
							      "m" (O)   //6
								      :
									      // register clobber list
									      "rax", "rcx", "rdx", "rsi", "rdi", "memory"
										      );
				      }
				      if (n0_rem)
				      {
					      R x0 = I[i0 * is0 + i1 * is1];
					      R x1 = I[i0 * is0 + i1 * is1 + 1];
					      O[i0 * os0 + i1 * os1] = x0;
					      O[i0 * os0 + i1 * os1 + 1] = x1;
				      }
			      }
			      break;

			      case 2://only os0 is 2. 256-bit contiguous write possible
			      n0 = n0 - n0_rem;
			      for (i1 = 0; i1 < n1; ++i1) {
				      for (i0 = 0; i0 < n0; i0+=2) {
					      in1_128 = _mm_loadu_pd((double const *)&I[i0 * is0 + i1 * is1]);
					      in2_128 = _mm_loadu_pd((double const *)&I[(i0+1) * is0 + i1 * is1]);
					      in1 = _mm256_castpd128_pd256(in1_128);
					      //in2 = _mm256_castpd128_pd256(in2_128);
					      //out1 = _mm256_permute2f128_pd(in1, in2, 0x20);
					      out1 = _mm256_insertf128_pd(in1, in2_128, 1);
					      _mm256_storeu_pd((double *)&O[i0 * os0 + i1 * os1], out1);
				      }
				      if (n0_rem)
				      {
					      R x0 = I[i0 * is0 + i1 * is1];
					      R x1 = I[i0 * is0 + i1 * is1 + 1];
					      O[i0 * os0 + i1 * os1] = x0;
					      O[i0 * os0 + i1 * os1 + 1] = x1;
				      }
			      }
			      break;

			      case 6://os0=2 and is1=2. Both 256-bit read and 256-bit write possible
			      n0 = n0 - n0_rem;
			      n1 = n1 - n1_rem;
			      for (i1 = 0; i1 < n1; i1+=2)
			      {
				      for (i0 = 0; i0 < n0; i0+=2) {
					      in1 = _mm256_loadu_pd((double const *)&I[i0 * is0 + i1 * is1]);
					      in2 = _mm256_loadu_pd((double const *)&I[(i0+1) * is0 + i1 * is1]);

					      //out1 = _mm256_shuffle_pd(in1, in2, 0x33);
					      //out2 = _mm256_shuffle_pd(in1, in2, 0x11);
					      out1 = _mm256_permute2f128_pd(in1, in2, 0x20);
					      out2 = _mm256_permute2f128_pd(in1, in2, 0x31);
					      _mm256_storeu_pd((double *)&O[i0 * os0 + i1 * os1], out1);
					      _mm256_storeu_pd((double *)&O[i0 * os0 + (i1+1) * os1], out2);
				      }
				      if (n0_rem)
				      {
					      R x0 = I[i0 * is0 + i1 * is1];
					      R x1 = I[i0 * is0 + i1 * is1 + 1];
					      R x2 = I[i0 * is0 + (i1+1) * is1];
					      R x3 = I[i0 * is0 + (i1+1) * is1 + 1];
					      O[i0 * os0 + i1 * os1] = x0;
					      O[i0 * os0 + i1 * os1 + 1] = x1;
					      O[i0 * os0 + (i1+1) * os1] = x2;
					      O[i0 * os0 + (i1+1) * os1 + 1] = x3;
				      }
			      }
			      n0 += n0_rem;
			      if (n1_rem)
			      {
				      for (i0 = 0; i0 < n0; ++i0) {
					      R x0 = I[i0 * is0 + i1 * is1];
					      R x1 = I[i0 * is0 + i1 * is1 + 1];
					      O[i0 * os0 + i1 * os1] = x0;
					      O[i0 * os0 + i1 * os1 + 1] = x1;
				      }
			      }
			      break;

			      case 8://only os1 is 2. 256-bit contiguous write possible
			      n1 = n1 - n1_rem;
			      _mm256_zeroupper();
			      for (i1 = 0; i1 < n1; i1+=2) {
				      t1 = i1 * is1;
				      t2 = i1 * os1;
				      for (i0 = 0; i0 < n0; ++i0) {
					      __asm__ volatile (
							      //"vzeroall \n\t"
							      "mov %0, %%rcx \n\t"
							      "mov %1, %%rax \n\t"
							      "mul %%rcx \n\t"
							      "mov %2, %%rdx \n\t"
							      "mov %5, %%rsi \n\t"
							      "add %3, %%rsi\n\t"
							      "mov %6, %%rdi\n\t"
							      "add %4, %%rdi\n\t"
							      "vmovupd (%%rsi, %%rax), %%xmm0\n\t"
							      "add %7, %%rax\n\t"
							      "vmovupd (%%rsi, %%rax), %%xmm1\n\t"
							      "mov %%rcx, %%rax \n\t"
							      "mul %%rdx \n\t"
							      "vinsertf128 $1, %%xmm0, %%ymm1, %%ymm2\n\t"
							      "vmovupd %%ymm2, (%%rdi, %%rax)\n\t"
							      :// output operands (none)
							      :// input operands
							      "m" (i0), //0
							      "m" (is0),//1
							      "m" (os0),//2
							      "m" (t1), //3
							      "m" (t2), //4
							      "m" (I),  //5
							      "m" (O),  //6
							      "m" (is1) //7
								      :
									      // register clobber list
									      "rax", "rcx", "rdx", "rsi", "rdi", "memory"
										      );
				      }
			      }
			      if (n1_rem)
			      {
				      for (i0 = 0; i0 < n0; ++i0) {
					      R x0 = I[i0 * is0 + i1 * is1];
					      R x1 = I[i0 * is0 + i1 * is1 + 1];
					      O[i0 * os0 + i1 * os1] = x0;
					      O[i0 * os0 + i1 * os1 + 1] = x1;
				      }
			      }
			      break;
			      
			      case 9://is0=2 and os1=2. Both 256-bit read and 256-bit write possible
			      n0 = n0 - n0_rem;
			      n1 = n1 - n1_rem;
			      for (i1 = 0; i1 < n1; i1+=2)
			      {
				      for (i0 = 0; i0 < n0; i0+=2) {
					      in1 = _mm256_loadu_pd((double const *)&I[i0 * is0 + i1 * is1]);
					      in2 = _mm256_loadu_pd((double const *)&I[i0 * is0 + (i1+1) * is1]);

					      //out1 = _mm256_shuffle_pd(in1, in2, 0x33);
					      //out2 = _mm256_shuffle_pd(in1, in2, 0x11);
					      out1 = _mm256_permute2f128_pd(in1, in2, 0x20);
					      out2 = _mm256_permute2f128_pd(in1, in2, 0x31);
					      _mm256_storeu_pd((double *)&O[i0 * os0 + i1 * os1], out1);
					      _mm256_storeu_pd((double *)&O[(i0+1) * os0 + i1 * os1], out2);
				      }
				      if (n0_rem)
				      {
					      R x0 = I[i0 * is0 + i1 * is1];
					      R x1 = I[i0 * is0 + i1 * is1 + 1];
					      R x2 = I[i0 * is0 + (i1+1) * is1];
					      R x3 = I[i0 * is0 + (i1+1) * is1 + 1];
					      O[i0 * os0 + i1 * os1] = x0;
					      O[i0 * os0 + i1 * os1 + 1] = x1;
					      O[i0 * os0 + (i1+1) * os1] = x2;
					      O[i0 * os0 + (i1+1) * os1 + 1] = x3;
				      }
			      }
			      if (n1_rem)
			      {
				      n0 += n0_rem;
				      for (i0 = 0; i0 < n0; ++i0) {
					      R x0 = I[i0 * is0 + i1 * is1];
					      R x1 = I[i0 * is0 + i1 * is1 + 1];
					      O[i0 * os0 + i1 * os1] = x0;
					      O[i0 * os0 + i1 * os1 + 1] = x1;
				      }
			      }
			      break;

			      case 3://is0=2 and os0=2. Both 256-bit read and 256-bit write possible
			      case 7://is0=2 and os0=2. Also is1=2. Both 256-bit read and 256-bit write possible
			      case 11://is0=2 and os0=2. Also os1=2. Both 256-bit read and 256-bit write possible
			      case 15://is0=2 and os0=2. Also is1=2, os1=2. Both 256-bit read and 256-bit write possible
			      n0 = n0 - n0_rem;
			      for (i1 = 0; i1 < n1; ++i1)
			      {
				      for (i0 = 0; i0 < n0; i0+=2) {
					      in1 = _mm256_loadu_pd((double const *)&I[i0 * is0 + i1 * is1]);
					      _mm256_storeu_pd((double *)&O[i0 * os0 + i1 * os1], in1);
				      }
				      if (n0_rem)
				      {
					      R x0 = I[i0 * is0 + i1 * is1];
					      R x1 = I[i0 * is0 + i1 * is1 + 1];
					      O[i0 * os0 + i1 * os1] = x0;
					      O[i0 * os0 + i1 * os1 + 1] = x1;
				      }
			      }
			      break;
			      
			      default:
			      for (i1 = 0; i1 < n1; ++i1)
				      for (i0 = 0; i0 < n0; ++i0) {
					      R x0 = I[i0 * is0 + i1 * is1];
					      R x1 = I[i0 * is0 + i1 * is1 + 1];
					      O[i0 * os0 + i1 * os1] = x0;
					      O[i0 * os0 + i1 * os1 + 1] = x1;
				      }
			      break;
		      }
	      }
	      break;
	 default:
	      for (i1 = 0; i1 < n1; ++i1)
		   for (i0 = 0; i0 < n0; ++i0)
			for (v = 0; v < vl; ++v) {
			     R x0 = I[i0 * is0 + i1 * is1 + v];
			     O[i0 * os0 + i1 * os1 + v] = x0;
			}
	      break;
     }
}
#endif//AMD_OPT_IN_PLACE_1D_CPY2D_EXPERIMENTAL_ASM

#else //original cpy2d routine

void X(cpy2d)(R *I, R *O,
	      INT n0, INT is0, INT os0,
	      INT n1, INT is1, INT os1,
	      INT vl)
{
     INT i0, i1, v;

     switch (vl) {
	 case 1:
	      for (i1 = 0; i1 < n1; ++i1)
		   for (i0 = 0; i0 < n0; ++i0) {
			R x0 = I[i0 * is0 + i1 * is1];
			O[i0 * os0 + i1 * os1] = x0;
		   }
	      break;
	 case 2:
	      if (1
		  && (2 * sizeof(R) == sizeof(WIDE_TYPE))
		  && (sizeof(WIDE_TYPE) > sizeof(double))
		  && (((size_t)I) % sizeof(WIDE_TYPE) == 0)
		  && (((size_t)O) % sizeof(WIDE_TYPE) == 0)
		  && ((is0 & 1) == 0)
		  && ((is1 & 1) == 0)
		  && ((os0 & 1) == 0)
		  && ((os1 & 1) == 0)) {
		   /* copy R[2] as WIDE_TYPE if WIDE_TYPE is large
		      enough to hold R[2], and if the input is
		      properly aligned.  This is a win when R==double
		      and WIDE_TYPE is 128 bits. */
		   for (i1 = 0; i1 < n1; ++i1)
			for (i0 = 0; i0 < n0; ++i0) {
			     *(WIDE_TYPE *)&O[i0 * os0 + i1 * os1] =
				  *(WIDE_TYPE *)&I[i0 * is0 + i1 * is1];
			}
	      } else if (1
		  && (2 * sizeof(R) == sizeof(double))
		  && (((size_t)I) % sizeof(double) == 0)
		  && (((size_t)O) % sizeof(double) == 0)
		  && ((is0 & 1) == 0)
		  && ((is1 & 1) == 0)
		  && ((os0 & 1) == 0)
		  && ((os1 & 1) == 0)) {
		   /* copy R[2] as double if double is large enough to
		      hold R[2], and if the input is properly aligned.
		      This case applies when R==float */
		   for (i1 = 0; i1 < n1; ++i1)
			for (i0 = 0; i0 < n0; ++i0) {
			     *(double *)&O[i0 * os0 + i1 * os1] =
				  *(double *)&I[i0 * is0 + i1 * is1];
			}
	      } else {
		   for (i1 = 0; i1 < n1; ++i1)
			for (i0 = 0; i0 < n0; ++i0) {
			     R x0 = I[i0 * is0 + i1 * is1];
			     R x1 = I[i0 * is0 + i1 * is1 + 1];
			     O[i0 * os0 + i1 * os1] = x0;
 			     O[i0 * os0 + i1 * os1 + 1] = x1;
			}
	      }
	      break;
	 default:
	      for (i1 = 0; i1 < n1; ++i1)
		   for (i0 = 0; i0 < n0; ++i0)
			for (v = 0; v < vl; ++v) {
			     R x0 = I[i0 * is0 + i1 * is1 + v];
			     O[i0 * os0 + i1 * os1 + v] = x0;
			}
	      break;
     }
}
#endif//AMD_OPT_ALL

/* like cpy2d, but read input contiguously if possible */
void X(cpy2d_ci)(R *I, R *O,
		 INT n0, INT is0, INT os0,
		 INT n1, INT is1, INT os1,
		 INT vl)
{
     if (IABS(is0) < IABS(is1))	/* inner loop is for n0 */
	  X(cpy2d) (I, O, n0, is0, os0, n1, is1, os1, vl);
     else
	  X(cpy2d) (I, O, n1, is1, os1, n0, is0, os0, vl);
}

/* like cpy2d, but write output contiguously if possible */
void X(cpy2d_co)(R *I, R *O,
		 INT n0, INT is0, INT os0,
		 INT n1, INT is1, INT os1,
		 INT vl)
{
     if (IABS(os0) < IABS(os1))	/* inner loop is for n0 */
	  X(cpy2d) (I, O, n0, is0, os0, n1, is1, os1, vl);
     else
	  X(cpy2d) (I, O, n1, is1, os1, n0, is0, os0, vl);
}


/* tiled copy routines */
struct cpy2d_closure {
     R *I, *O;
     INT is0, os0, is1, os1, vl;
     R *buf;
};

static void dotile(INT n0l, INT n0u, INT n1l, INT n1u, void *args)
{
     struct cpy2d_closure *k = (struct cpy2d_closure *)args;
     X(cpy2d)(k->I + n0l * k->is0 + n1l * k->is1,
	      k->O + n0l * k->os0 + n1l * k->os1,
	      n0u - n0l, k->is0, k->os0,
	      n1u - n1l, k->is1, k->os1,
	      k->vl);
}

static void dotile_buf(INT n0l, INT n0u, INT n1l, INT n1u, void *args)
{
     struct cpy2d_closure *k = (struct cpy2d_closure *)args;

     /* copy from I to buf */
     X(cpy2d_ci)(k->I + n0l * k->is0 + n1l * k->is1,
		 k->buf,
		 n0u - n0l, k->is0, k->vl,
		 n1u - n1l, k->is1, k->vl * (n0u - n0l),
		 k->vl);

     /* copy from buf to O */
     X(cpy2d_co)(k->buf,
		 k->O + n0l * k->os0 + n1l * k->os1,
		 n0u - n0l, k->vl, k->os0,
		 n1u - n1l, k->vl * (n0u - n0l), k->os1,
		 k->vl);
}


void X(cpy2d_tiled)(R *I, R *O,
		    INT n0, INT is0, INT os0,
		    INT n1, INT is1, INT os1, INT vl)
{
     INT tilesz = X(compute_tilesz)(vl,
				    1 /* input array */
				    + 1 /* ouput array */);
     struct cpy2d_closure k;
     k.I = I;
     k.O = O;
     k.is0 = is0;
     k.os0 = os0;
     k.is1 = is1;
     k.os1 = os1;
     k.vl = vl;
     k.buf = 0; /* unused */
     X(tile2d)(0, n0, 0, n1, tilesz, dotile, &k);
}

void X(cpy2d_tiledbuf)(R *I, R *O,
		       INT n0, INT is0, INT os0,
		       INT n1, INT is1, INT os1, INT vl)
{
     R buf[CACHESIZE / (2 * sizeof(R))];
     /* input and buffer in cache, or
	output and buffer in cache */
     INT tilesz = X(compute_tilesz)(vl, 2);
     struct cpy2d_closure k;
     k.I = I;
     k.O = O;
     k.is0 = is0;
     k.os0 = os0;
     k.is1 = is1;
     k.os1 = os1;
     k.vl = vl;
     k.buf = buf;
     A(tilesz * tilesz * vl * sizeof(R) <= sizeof(buf));
     X(tile2d)(0, n0, 0, n1, tilesz, dotile_buf, &k);
}
