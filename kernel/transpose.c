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

#include "kernel/ifftw.h"

#if 0 //#ifdef AMD_FMV_AUTO //Let the manual FMV option be enabled
/* in place square transposition, iterative */
__attribute__((target_clones(TARGET_STRINGS)))
void X(transpose)(R *I, INT n, INT s0, INT s1, INT vl)
{
     INT i0, i1, v;

     switch (vl) {
	 case 1:
          for (i1 = 1; i1 < n; ++i1) {
           for (i0 = 0; i0 < i1; ++i0) {
            R x0 = I[i1 * s0 + i0 * s1];
            R y0 = I[i1 * s1 + i0 * s0];
            I[i1 * s1 + i0 * s0] = x0;
            I[i1 * s0 + i0 * s1] = y0;
           }
          }
          break;
	 case 2:
          for (i1 = 1; i1 < n; ++i1) {
           for (i0 = 0; i0 < i1; ++i0) {
            R x0 = I[i1 * s0 + i0 * s1];
            R x1 = I[i1 * s0 + i0 * s1 + 1];
            R y0 = I[i1 * s1 + i0 * s0];
            R y1 = I[i1 * s1 + i0 * s0 + 1];
            I[i1 * s1 + i0 * s0] = x0;
            I[i1 * s1 + i0 * s0 + 1] = x1;
            I[i1 * s0 + i0 * s1] = y0;
            I[i1 * s0 + i0 * s1 + 1] = y1;
           }
          }
          break;
	 default:
          for (i1 = 1; i1 < n; ++i1) {
           for (i0 = 0; i0 < i1; ++i0) {
            for (v = 0; v < vl; ++v) {
                R x0 = I[i1 * s0 + i0 * s1 + v];
                R y0 = I[i1 * s1 + i0 * s0 + v];
                I[i1 * s1 + i0 * s0 + v] = x0;
                I[i1 * s0 + i0 * s1 + v] = y0;
            }
           }
          }
          break;
     }
}
#else
#if defined(AMD_OPT_IN_PLACE_SQU_TRANS) && (!defined(FFTW_LDOUBLE) && !defined(FFTW_QUAD))
#include "immintrin.h"
#endif

#if defined(AMD_OPT_AUTO_TUNED_RASTER_TILED_TRANS_METHOD)
unsigned int L1D_blk_size;// = CACHESIZE;
#endif

/* in place square transposition, iterative */

#ifdef AMD_FMV_MANUAL
//variant of transpose routine for AVX
__attribute__ ((target ("avx")))
void X(transpose_avx)(R *I, INT n, INT s0, INT s1, INT vl)
{
     INT i0, i1, v;
#if defined(AMD_OPT_IN_PLACE_SQU_TRANS) && (!defined(FFTW_LDOUBLE) && !defined(FFTW_QUAD))
	 int vl_8factor = vl - (vl & 0xF);
#ifdef FFTW_SINGLE
 	 __m256 in1, in2, in3, in4, in5, in6, in7, in8;
#else
	 __m256d in1, in2, in3, in4, in5, in6, in7, in8;
#endif
#endif

     switch (vl) {
	 case 1:
	      for (i1 = 1; i1 < n; ++i1) {
		   for (i0 = 0; i0 < i1; ++i0) {
			R x0 = I[i1 * s0 + i0 * s1];
			R y0 = I[i1 * s1 + i0 * s0];
			I[i1 * s1 + i0 * s0] = x0;
			I[i1 * s0 + i0 * s1] = y0;
		   }
	      }
	      break;
	 case 2:
	      for (i1 = 1; i1 < n; ++i1) {
		   for (i0 = 0; i0 < i1; ++i0) {
			R x0 = I[i1 * s0 + i0 * s1];
			R x1 = I[i1 * s0 + i0 * s1 + 1];
			R y0 = I[i1 * s1 + i0 * s0];
			R y1 = I[i1 * s1 + i0 * s0 + 1];
			I[i1 * s1 + i0 * s0] = x0;
			I[i1 * s1 + i0 * s0 + 1] = x1;
			I[i1 * s0 + i0 * s1] = y0;
			I[i1 * s0 + i0 * s1 + 1] = y1;
		   }
	      }
	      break;
	 default:
#if defined(AMD_OPT_IN_PLACE_SQU_TRANS) && (!defined(FFTW_LDOUBLE) && !defined(FFTW_QUAD))
#ifdef FFTW_SINGLE
	      for (i1 = 1; i1 < n; ++i1) {
		      for (i0 = 0; i0 < i1; ++i0) {
			      for (v = 0; v < vl_8factor; v+=16) {
				      in1 = _mm256_loadu_ps((float const *)&I[i1 * s0 + i0 * s1 + v + 0]);
				      in2 = _mm256_loadu_ps((float const *)&I[i1 * s1 + i0 * s0 + v + 0]);
				      in3 = _mm256_loadu_ps((float const *)&I[i1 * s0 + i0 * s1 + v + 8]);
				      in4 = _mm256_loadu_ps((float const *)&I[i1 * s1 + i0 * s0 + v + 8]);
				      _mm256_storeu_ps((double *)&I[i1 * s1 + i0 * s0 + v + 0], in1);
				      _mm256_storeu_ps((double *)&I[i1 * s0 + i0 * s1 + v + 0], in2);
				      _mm256_storeu_ps((double *)&I[i1 * s1 + i0 * s0 + v + 8], in3);
				      _mm256_storeu_ps((double *)&I[i1 * s0 + i0 * s1 + v + 8], in4);
			      }
			      for (; v < vl; ++v) {
				      R x0 = I[i1 * s0 + i0 * s1 + v];
				      R y0 = I[i1 * s1 + i0 * s0 + v];
				      I[i1 * s1 + i0 * s0 + v] = x0;
				      I[i1 * s0 + i0 * s1 + v] = y0;
			      }
		      }
	      }
#else
	      for (i1 = 1; i1 < n; ++i1) {
		      for (i0 = 0; i0 < i1; ++i0) {
			      for (v = 0; v < vl_8factor; v+=16) {
				      in1 = _mm256_loadu_pd((double const *)&I[i1 * s0 + i0 * s1 + v + 0]);
				      in2 = _mm256_loadu_pd((double const *)&I[i1 * s1 + i0 * s0 + v + 0]);
				      in3 = _mm256_loadu_pd((double const *)&I[i1 * s0 + i0 * s1 + v + 4]);
				      in4 = _mm256_loadu_pd((double const *)&I[i1 * s1 + i0 * s0 + v + 4]);
				      in5 = _mm256_loadu_pd((double const *)&I[i1 * s0 + i0 * s1 + v + 8]);
				      in6 = _mm256_loadu_pd((double const *)&I[i1 * s1 + i0 * s0 + v + 8]);
				      in7 = _mm256_loadu_pd((double const *)&I[i1 * s0 + i0 * s1 + v + 12]);
				      in8 = _mm256_loadu_pd((double const *)&I[i1 * s1 + i0 * s0 + v + 12]);
				      _mm256_storeu_pd((double *)&I[i1 * s1 + i0 * s0 + v + 0], in1);
				      _mm256_storeu_pd((double *)&I[i1 * s0 + i0 * s1 + v + 0], in2);
				      _mm256_storeu_pd((double *)&I[i1 * s1 + i0 * s0 + v + 4], in3);
				      _mm256_storeu_pd((double *)&I[i1 * s0 + i0 * s1 + v + 4], in4);
				      _mm256_storeu_pd((double *)&I[i1 * s1 + i0 * s0 + v + 8], in5);
				      _mm256_storeu_pd((double *)&I[i1 * s0 + i0 * s1 + v + 8], in6);
				      _mm256_storeu_pd((double *)&I[i1 * s1 + i0 * s0 + v + 12], in7);
				      _mm256_storeu_pd((double *)&I[i1 * s0 + i0 * s1 + v + 12], in8);
			      }
			      for (; v < vl; ++v) {
				      R x0 = I[i1 * s0 + i0 * s1 + v];
				      R y0 = I[i1 * s1 + i0 * s0 + v];
				      I[i1 * s1 + i0 * s0 + v] = x0;
				      I[i1 * s0 + i0 * s1 + v] = y0;
			      }
		      }
	      }
#endif
#else
	      for (i1 = 1; i1 < n; ++i1) {
		   for (i0 = 0; i0 < i1; ++i0) {
			for (v = 0; v < vl; ++v) {
			     R x0 = I[i1 * s0 + i0 * s1 + v];
			     R y0 = I[i1 * s1 + i0 * s0 + v];
			     I[i1 * s1 + i0 * s0 + v] = x0;
			     I[i1 * s0 + i0 * s1 + v] = y0;
			}
		   }
	      }
#endif
	      break;
     }
}

//variant of transpose routine for AVX512
#ifdef HAVE_AVX512
__attribute__ ((target ("avx512f")))
void X(transpose_avx512)(R *I, INT n, INT s0, INT s1, INT vl)
{
	INT i0, i1, v;
#if defined(AMD_OPT_IN_PLACE_SQU_TRANS) && (!defined(FFTW_LDOUBLE) && !defined(FFTW_QUAD))
		int vl_8factor = vl - (vl & 0xF);
#ifdef FFTW_SINGLE
		__m512 in1, in2, in3, in4;
#else
		__m512d in1, in2, in3, in4;
#endif
#endif
		switch (vl) {
			case 1:
				for (i1 = 1; i1 < n; ++i1) {
					for (i0 = 0; i0 < i1; ++i0) {
						R x0 = I[i1 * s0 + i0 * s1];
						R y0 = I[i1 * s1 + i0 * s0];
						I[i1 * s1 + i0 * s0] = x0;
						I[i1 * s0 + i0 * s1] = y0;
					}
				}
				break;
			case 2:
				for (i1 = 1; i1 < n; ++i1) {
					for (i0 = 0; i0 < i1; ++i0) {
						R x0 = I[i1 * s0 + i0 * s1];
						R x1 = I[i1 * s0 + i0 * s1 + 1];
						R y0 = I[i1 * s1 + i0 * s0];
						R y1 = I[i1 * s1 + i0 * s0 + 1];
						I[i1 * s1 + i0 * s0] = x0;
						I[i1 * s1 + i0 * s0 + 1] = x1;
						I[i1 * s0 + i0 * s1] = y0;
						I[i1 * s0 + i0 * s1 + 1] = y1;
					}
				}
				break;
			default:
#if defined(AMD_OPT_IN_PLACE_SQU_TRANS) && (!defined(FFTW_LDOUBLE) && !defined(FFTW_QUAD))
#ifdef FFTW_SINGLE
				for (i1 = 1; i1 < n; ++i1) {
					for (i0 = 0; i0 < i1; ++i0) {
						for (v = 0; v < vl_8factor; v+=16) {
							in1 = _mm512_loadu_ps((float const *)&I[i1 * s0 + i0 * s1 + v + 0]);
							in2 = _mm512_loadu_ps((float const *)&I[i1 * s1 + i0 * s0 + v + 0]);
							_mm512_storeu_ps((double *)&I[i1 * s1 + i0 * s0 + v + 0], in1);
							_mm512_storeu_ps((double *)&I[i1 * s0 + i0 * s1 + v + 0], in2);
						}
						for (; v < vl; ++v) {
							R x0 = I[i1 * s0 + i0 * s1 + v];
							R y0 = I[i1 * s1 + i0 * s0 + v];
							I[i1 * s1 + i0 * s0 + v] = x0;
							I[i1 * s0 + i0 * s1 + v] = y0;
						}
					}
				}
#else
				for (i1 = 1; i1 < n; ++i1) {
					for (i0 = 0; i0 < i1; ++i0) {
						for (v = 0; v < vl_8factor; v+=16) {
							in1 = _mm512_loadu_pd((double const *)&I[i1 * s0 + i0 * s1 + v + 0]);
							in2 = _mm512_loadu_pd((double const *)&I[i1 * s1 + i0 * s0 + v + 0]);
							in3 = _mm512_loadu_pd((double const *)&I[i1 * s0 + i0 * s1 + v + 8]);
							in4 = _mm512_loadu_pd((double const *)&I[i1 * s1 + i0 * s0 + v + 8]);
							_mm512_storeu_pd((double *)&I[i1 * s1 + i0 * s0 + v + 0], in1);
							_mm512_storeu_pd((double *)&I[i1 * s0 + i0 * s1 + v + 0], in2);
							_mm512_storeu_pd((double *)&I[i1 * s1 + i0 * s0 + v + 8], in3);
							_mm512_storeu_pd((double *)&I[i1 * s0 + i0 * s1 + v + 8], in4);
						}
						for (; v < vl; ++v) {
							R x0 = I[i1 * s0 + i0 * s1 + v];
							R y0 = I[i1 * s1 + i0 * s0 + v];
							I[i1 * s1 + i0 * s0 + v] = x0;
							I[i1 * s0 + i0 * s1 + v] = y0;
						}
					}
				}
#endif
#else
				for (i1 = 1; i1 < n; ++i1) {
					for (i0 = 0; i0 < i1; ++i0) {
						for (v = 0; v < vl; ++v) {
							R x0 = I[i1 * s0 + i0 * s1 + v];
							R y0 = I[i1 * s1 + i0 * s0 + v];
							I[i1 * s1 + i0 * s0 + v] = x0;
							I[i1 * s0 + i0 * s1 + v] = y0;
						}
					}
				}
#endif
				break;
		}
}
#endif
__attribute__ ((target ("sse2")))
void X(transpose_c)(R *I, INT n, INT s0, INT s1, INT vl)
{
     INT i0, i1, v;

     switch (vl) {
	 case 1:
          for (i1 = 1; i1 < n; ++i1) {
           for (i0 = 0; i0 < i1; ++i0) {
            R x0 = I[i1 * s0 + i0 * s1];
            R y0 = I[i1 * s1 + i0 * s0];
            I[i1 * s1 + i0 * s0] = x0;
            I[i1 * s0 + i0 * s1] = y0;
           }
          }
          break;
	 case 2:
          for (i1 = 1; i1 < n; ++i1) {
           for (i0 = 0; i0 < i1; ++i0) {
            R x0 = I[i1 * s0 + i0 * s1];
            R x1 = I[i1 * s0 + i0 * s1 + 1];
            R y0 = I[i1 * s1 + i0 * s0];
            R y1 = I[i1 * s1 + i0 * s0 + 1];
            I[i1 * s1 + i0 * s0] = x0;
            I[i1 * s1 + i0 * s0 + 1] = x1;
            I[i1 * s0 + i0 * s1] = y0;
            I[i1 * s0 + i0 * s1 + 1] = y1;
           }
          }
          break;
	 default:
          for (i1 = 1; i1 < n; ++i1) {
           for (i0 = 0; i0 < i1; ++i0) {
            for (v = 0; v < vl; ++v) {
                R x0 = I[i1 * s0 + i0 * s1 + v];
                R y0 = I[i1 * s1 + i0 * s0 + v];
                I[i1 * s1 + i0 * s0 + v] = x0;
                I[i1 * s0 + i0 * s1 + v] = y0;
            }
           }
          }
          break;
     }
}


int is_avx512, is_avx;

void X(transpose_internal)(R* I, INT n, INT s0, INT s1, INT vl) __attribute__((ifunc("fmv_resolver_transpose")));
void X(transpose)(R* I, INT n, INT s0, INT s1, INT vl)
{
#if defined(HAVE_AVX512)
	is_avx512 = X(have_simd_avx512);
#else
	is_avx512 = 0;
#endif
#if defined(HAVE_AVX)
	is_avx = X(have_simd_avx);
#else
	is_avx = 0;
#endif
	X(transpose_internal)(I, n, s0, s1, vl);
}

static void* fmv_resolver_transpose(void)
{
#if defined(HAVE_AVX512)
	if (is_avx512)
	{
		return X(transpose_avx512);
	}
	else
#endif
	{
#if defined(HAVE_AVX)
		if (is_avx)
		{
			return X(transpose_avx);
		}
		else
#endif
		{
			return X(transpose_c);
		}
	}

}
#else
void X(transpose)(R *I, INT n, INT s0, INT s1, INT vl)
{
     INT i0, i1, v;
#if defined(AMD_OPT_IN_PLACE_SQU_TRANS) && (!defined(FFTW_LDOUBLE) && !defined(FFTW_QUAD))
	 int vl_8factor = vl - (vl & 0xF);
#ifdef FFTW_SINGLE
 	 __m256 in1, in2, in3, in4, in5, in6, in7, in8;
#else
	 __m256d in1, in2, in3, in4, in5, in6, in7, in8;
#endif
#endif

     switch (vl) {
	 case 1:
	      for (i1 = 1; i1 < n; ++i1) {
		   for (i0 = 0; i0 < i1; ++i0) {
			R x0 = I[i1 * s0 + i0 * s1];
			R y0 = I[i1 * s1 + i0 * s0];
			I[i1 * s1 + i0 * s0] = x0;
			I[i1 * s0 + i0 * s1] = y0;
		   }
	      }
	      break;
	 case 2:
	      for (i1 = 1; i1 < n; ++i1) {
		   for (i0 = 0; i0 < i1; ++i0) {
			R x0 = I[i1 * s0 + i0 * s1];
			R x1 = I[i1 * s0 + i0 * s1 + 1];
			R y0 = I[i1 * s1 + i0 * s0];
			R y1 = I[i1 * s1 + i0 * s0 + 1];
			I[i1 * s1 + i0 * s0] = x0;
			I[i1 * s1 + i0 * s0 + 1] = x1;
			I[i1 * s0 + i0 * s1] = y0;
			I[i1 * s0 + i0 * s1 + 1] = y1;
		   }
	      }
	      break;
	 default:
#if defined(AMD_OPT_IN_PLACE_SQU_TRANS) && (!defined(FFTW_LDOUBLE) && !defined(FFTW_QUAD))
#ifdef FFTW_SINGLE
	      for (i1 = 1; i1 < n; ++i1) {
		      for (i0 = 0; i0 < i1; ++i0) {
			      for (v = 0; v < vl_8factor; v+=16) {
				      in1 = _mm256_loadu_ps((float const *)&I[i1 * s0 + i0 * s1 + v + 0]);
				      in2 = _mm256_loadu_ps((float const *)&I[i1 * s1 + i0 * s0 + v + 0]);
				      in3 = _mm256_loadu_ps((float const *)&I[i1 * s0 + i0 * s1 + v + 8]);
				      in4 = _mm256_loadu_ps((float const *)&I[i1 * s1 + i0 * s0 + v + 8]);
				      _mm256_storeu_ps((double *)&I[i1 * s1 + i0 * s0 + v + 0], in1);
				      _mm256_storeu_ps((double *)&I[i1 * s0 + i0 * s1 + v + 0], in2);
				      _mm256_storeu_ps((double *)&I[i1 * s1 + i0 * s0 + v + 8], in3);
				      _mm256_storeu_ps((double *)&I[i1 * s0 + i0 * s1 + v + 8], in4);
			      }
			      for (; v < vl; ++v) {
				      R x0 = I[i1 * s0 + i0 * s1 + v];
				      R y0 = I[i1 * s1 + i0 * s0 + v];
				      I[i1 * s1 + i0 * s0 + v] = x0;
				      I[i1 * s0 + i0 * s1 + v] = y0;
			      }
		      }
	      }
#else
	      for (i1 = 1; i1 < n; ++i1) {
		      for (i0 = 0; i0 < i1; ++i0) {
			      for (v = 0; v < vl_8factor; v+=16) {
				      in1 = _mm256_loadu_pd((double const *)&I[i1 * s0 + i0 * s1 + v + 0]);
				      in2 = _mm256_loadu_pd((double const *)&I[i1 * s1 + i0 * s0 + v + 0]);
				      in3 = _mm256_loadu_pd((double const *)&I[i1 * s0 + i0 * s1 + v + 4]);
				      in4 = _mm256_loadu_pd((double const *)&I[i1 * s1 + i0 * s0 + v + 4]);
				      in5 = _mm256_loadu_pd((double const *)&I[i1 * s0 + i0 * s1 + v + 8]);
				      in6 = _mm256_loadu_pd((double const *)&I[i1 * s1 + i0 * s0 + v + 8]);
				      in7 = _mm256_loadu_pd((double const *)&I[i1 * s0 + i0 * s1 + v + 12]);
				      in8 = _mm256_loadu_pd((double const *)&I[i1 * s1 + i0 * s0 + v + 12]);
				      _mm256_storeu_pd((double *)&I[i1 * s1 + i0 * s0 + v + 0], in1);
				      _mm256_storeu_pd((double *)&I[i1 * s0 + i0 * s1 + v + 0], in2);
				      _mm256_storeu_pd((double *)&I[i1 * s1 + i0 * s0 + v + 4], in3);
				      _mm256_storeu_pd((double *)&I[i1 * s0 + i0 * s1 + v + 4], in4);
				      _mm256_storeu_pd((double *)&I[i1 * s1 + i0 * s0 + v + 8], in5);
				      _mm256_storeu_pd((double *)&I[i1 * s0 + i0 * s1 + v + 8], in6);
				      _mm256_storeu_pd((double *)&I[i1 * s1 + i0 * s0 + v + 12], in7);
				      _mm256_storeu_pd((double *)&I[i1 * s0 + i0 * s1 + v + 12], in8);
			      }
			      for (; v < vl; ++v) {
				      R x0 = I[i1 * s0 + i0 * s1 + v];
				      R y0 = I[i1 * s1 + i0 * s0 + v];
				      I[i1 * s1 + i0 * s0 + v] = x0;
				      I[i1 * s0 + i0 * s1 + v] = y0;
			      }
		      }
	      }
#endif
#else
	      for (i1 = 1; i1 < n; ++i1) {
		   for (i0 = 0; i0 < i1; ++i0) {
			for (v = 0; v < vl; ++v) {
			     R x0 = I[i1 * s0 + i0 * s1 + v];
			     R y0 = I[i1 * s1 + i0 * s0 + v];
			     I[i1 * s1 + i0 * s0 + v] = x0;
			     I[i1 * s0 + i0 * s1 + v] = y0;
			}
		   }
	      }
#endif
	      break;
     }
}
#endif
#endif

struct transpose_closure {
     R *I;
     INT s0, s1, vl, tilesz;
     R *buf0, *buf1; 
};

#ifdef AMD_FMV_AUTO
__attribute__((target_clones(TARGET_STRINGS)))
#endif
static void dotile(INT n0l, INT n0u, INT n1l, INT n1u, void *args)
{
     struct transpose_closure *k = (struct transpose_closure *)args;
     R *I = k->I;
     INT s0 = k->s0, s1 = k->s1, vl = k->vl;
     INT i0, i1, v;

     switch (vl) {
	 case 1:
	      for (i1 = n1l; i1 < n1u; ++i1) {
		   for (i0 = n0l; i0 < n0u; ++i0) {
			R x0 = I[i1 * s0 + i0 * s1];
			R y0 = I[i1 * s1 + i0 * s0];
			I[i1 * s1 + i0 * s0] = x0;
			I[i1 * s0 + i0 * s1] = y0;
		   }
	      }
	      break;
	 case 2:
	      for (i1 = n1l; i1 < n1u; ++i1) {
		   for (i0 = n0l; i0 < n0u; ++i0) {
			R x0 = I[i1 * s0 + i0 * s1];
			R x1 = I[i1 * s0 + i0 * s1 + 1];
			R y0 = I[i1 * s1 + i0 * s0];
			R y1 = I[i1 * s1 + i0 * s0 + 1];
			I[i1 * s1 + i0 * s0] = x0;
			I[i1 * s1 + i0 * s0 + 1] = x1;
			I[i1 * s0 + i0 * s1] = y0;
			I[i1 * s0 + i0 * s1 + 1] = y1;
		   }
	      }
	      break;
	 default:
	      for (i1 = n1l; i1 < n1u; ++i1) {
		   for (i0 = n0l; i0 < n0u; ++i0) {
			for (v = 0; v < vl; ++v) {
			     R x0 = I[i1 * s0 + i0 * s1 + v];
			     R y0 = I[i1 * s1 + i0 * s0 + v];
			     I[i1 * s1 + i0 * s0 + v] = x0;
			     I[i1 * s0 + i0 * s1 + v] = y0;
			}
		   }
	      }
     }
}

static void dotile_buf(INT n0l, INT n0u, INT n1l, INT n1u, void *args)
{
     struct transpose_closure *k = (struct transpose_closure *)args;
     X(cpy2d_ci)(k->I + n0l * k->s0 + n1l * k->s1,
		 k->buf0,
		 n0u - n0l, k->s0, k->vl,
		 n1u - n1l, k->s1, k->vl * (n0u - n0l),
		 k->vl);
     X(cpy2d_ci)(k->I + n0l * k->s1 + n1l * k->s0,
		 k->buf1,
		 n0u - n0l, k->s1, k->vl,
		 n1u - n1l, k->s0, k->vl * (n0u - n0l),
		 k->vl);
     X(cpy2d_co)(k->buf1,
		 k->I + n0l * k->s0 + n1l * k->s1,
		 n0u - n0l, k->vl, k->s0,
		 n1u - n1l, k->vl * (n0u - n0l), k->s1,
		 k->vl);
     X(cpy2d_co)(k->buf0,
		 k->I + n0l * k->s1 + n1l * k->s0,
		 n0u - n0l, k->vl, k->s1,
		 n1u - n1l, k->vl * (n0u - n0l), k->s0,
		 k->vl);
}

static void transpose_rec(R *I, INT n,
			  void (*f)(INT n0l, INT n0u, INT n1l, INT n1u,
				    void *args),
			  struct transpose_closure *k)
{
   tail:
     if (n > 1) {
	  INT n2 = n / 2;
	  k->I = I;
	  X(tile2d)(0, n2, n2, n, k->tilesz, f, k);
	  transpose_rec(I, n2, f, k);
	  I += n2 * (k->s0 + k->s1); n -= n2; goto tail;
     }
}

void X(transpose_tiled)(R *I, INT n, INT s0, INT s1, INT vl) 
{
     struct transpose_closure k;
     k.s0 = s0;
     k.s1 = s1;
     k.vl = vl;
     /* two blocks must be in cache, to be swapped */
     k.tilesz = X(compute_tilesz)(vl, 2);
     k.buf0 = k.buf1 = 0; /* unused */
     transpose_rec(I, n, dotile, &k);
}

#ifdef AMD_OPT_AUTO_TUNED_RASTER_TILED_TRANS_METHOD
void trans_autoTuned_tiled_Rorder(R *in, int n, struct transpose_closure *pk)
{
    int i, j, k, l;
    //static arrays of size equal to upper bound of L1D cache size, but it's utilized only as per current CPU's actual L1D size
#ifdef FFTW_SINGLE
    float buf1[2*BLK_SIZE*BLK_SIZE];
    float buf2[2*BLK_SIZE*BLK_SIZE];
#else
    double buf1[BLK_SIZE*BLK_SIZE];
    double buf2[BLK_SIZE*BLK_SIZE];
#endif
    int atc_blk_size = L1D_blk_size;//auto-tuned L1D cache block size
    //(vl >> 1) supports both double-precision and single-precision
    int num_ele = L1D_blk_size >> (pk->vl >> 1);//number of complex numbers is half for double-precision and one-fourth for single-precision.
    int m = pk->s0;
 
    //case when n is less than block size is handled by calling original fftw transpose function

    for (i = 0; i < n; i += num_ele)
    {
#ifdef FFTW_SINGLE
	float tmp1, tmp2;
#else
        double tmp1, tmp2;
#endif
	j = i;
	k = i<<1;
	//diagonal block: can be directly copied like ramModel
     	X(cpy2d_ci)(in + i * m + k,
		 buf1,
		 num_ele, pk->vl, pk->vl,
		 num_ele, pk->s0, atc_blk_size,
		 pk->vl);
     	X(cpy2d_co)(buf1,
		 in + i * m + k,
		 num_ele, atc_blk_size, pk->vl,
		 num_ele, pk->vl, pk->s0,
		 pk->vl);
	j += num_ele;

	//Next block in the raster scan order for current value of i
	for (; j < n; j += num_ele)
	{
		k = j << 1;
		l = i << 1;
		X(cpy2d_ci)(in + i * m + k,
				buf1,
				num_ele, pk->vl, pk->vl,
				num_ele, pk->s0, atc_blk_size,
				pk->vl);
		X(cpy2d_ci)(in + j * m + l,
				buf2,
				num_ele, pk->vl, pk->vl,
				num_ele, pk->s0, atc_blk_size,
				pk->vl);
		X(cpy2d_co)(buf2,
				in + i * m + k,
				num_ele, atc_blk_size, pk->vl,
				num_ele, pk->vl, pk->s0,
				pk->vl);
		X(cpy2d_co)(buf1,
				in + j * m + l,
				num_ele, atc_blk_size, pk->vl,
				num_ele, pk->vl, pk->s0,
				pk->vl);
	}
    }
}
#endif

void X(transpose_tiledbuf)(R *I, INT n, INT s0, INT s1, INT vl) 
{
     struct transpose_closure k;
#ifdef AMD_OPT_AUTO_TUNED_RASTER_TILED_TRANS_METHOD
     int blkSize = L1D_blk_size >> (vl >> 1);// (vl >> 1) supports both double and single precision
#endif
     /* Assume that the the rows of I conflict into the same cache
        lines, and therefore we don't need to reserve cache space for
        the input.  If the rows don't conflict, there is no reason
	to use tiledbuf at all.*/
     R buf0[CACHESIZE / (2 * sizeof(R))];
     R buf1[CACHESIZE / (2 * sizeof(R))];
     k.s0 = s0;
     k.s1 = s1;
     k.vl = vl;
     k.tilesz = X(compute_tilesz)(vl, 2);
     k.buf0 = buf0;
     k.buf1 = buf1;
     A(k.tilesz * k.tilesz * vl * sizeof(R) <= sizeof(buf0));
     A(k.tilesz * k.tilesz * vl * sizeof(R) <= sizeof(buf1));

#ifndef AMD_OPT_AUTO_TUNED_RASTER_TILED_TRANS_METHOD
     transpose_rec(I, n, dotile_buf, &k);
#else
     //Call original cache-oblivious transpose for cases:-
     //(i) when matrix size is smaller than block size
     //(ii) when matrix size is not multiple of block size
     //(iii) when vector length is 1. It is assumed that real input data will have vl=1 and enter here.
     if ((n < blkSize) || (n & (blkSize-1)) || (vl == 1))
     {
     	transpose_rec(I, n, dotile_buf, &k);
	return;
     }
     trans_autoTuned_tiled_Rorder(I, n, &k);
#endif
}

