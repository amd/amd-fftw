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

/* out of place copy routines for pairs of isomorphic 2D arrays */
#include "kernel/ifftw.h"

#ifdef AMD_OPT_ALL
#include "immintrin.h"
#endif

#if defined(AMD_OPT_ALL) && (!defined(FFTW_LDOUBLE) && !defined(FFTW_QUAD)) //AMD optimized routines

#ifdef FFTW_SINGLE//SINGLE PRECISION
#if defined(AMD_OPT_IN_PLACE_1D_CPY2D_STABLE_INTRIN)//SIMD optimized function
void X(cpy2d_pair)(R *I0, R *I1, R *O0, R *O1,
           INT n0, INT is0, INT os0,
           INT n1, INT is1, INT os1)
{
     INT i0, i1, v;
     R *I = I0, *O = O0;
     __m256 in1, in2, in3, in4;
     __m256 out1, out2, out3, out4;
     INT t0, t1, t2, t3;
     INT n0_rem, n1_rem;
     INT n0Many;
     INT n1Many;
     
     switch (((I1 - I0) == 1) & ((O1 - O0) == 1))
     {
        case 1:
        {
            n0_rem = n0&0x3, n1_rem = n1&0x3;
            n0Many = n0>3;
            n1Many = n1>3;
            t0 = (is0==2) & n0Many;
            t1 = (os0==2) & n0Many;
            t2 = (is1==2) & n1Many;
            t3 = (os1==2) & n1Many;
            
            switch(t0 | (t1 << 1) | (t2 << 2) | (t3 << 3))
            {
              case 6://os0=2 and is1=2. Both 256-bit read and 256-bit write possible
              n0 = n0 - n0_rem;
              n1 = n1 - n1_rem;
              for (i1 = 0; i1 < n1; i1+=4)
              {
                  for (i0 = 0; i0 < n0; i0+=4) 
                  {
                    in1 = _mm256_loadu_ps((float const *)&I[i0 * is0 + i1 * is1]);
                    in2 = _mm256_loadu_ps((float const *)&I[(i0+1) * is0 + i1 * is1]);
                    in3 = _mm256_loadu_ps((float const *)&I[(i0+2) * is0 + i1 * is1]);
                    in4 = _mm256_loadu_ps((float const *)&I[(i0+3) * is0 + i1 * is1]);
    
                    out2 = _mm256_shuffle_ps(in1, in2, 0x44);
                    out4 = _mm256_shuffle_ps(in3, in4, 0x44);
                    out1 = _mm256_permute2f128_ps(out2, out4, 0x20);
                    out3 = _mm256_permute2f128_ps(out2, out4, 0x31);
                    in1 = _mm256_shuffle_ps(in1, in2, 0xEE);
                    in3 = _mm256_shuffle_ps(in3, in4, 0xEE);
                    out2 = _mm256_permute2f128_ps(in1, in3, 0x20);
                    out4 = _mm256_permute2f128_ps(in1, in3, 0x31);
    
                    _mm256_storeu_ps((float *)&O[i0 * os0 + i1 * os1], out1);
                    _mm256_storeu_ps((float *)&O[i0 * os0 + (i1+1) * os1], out2);
                    _mm256_storeu_ps((float *)&O[i0 * os0 + (i1+2) * os1], out3);
                    _mm256_storeu_ps((float *)&O[i0 * os0 + (i1+3) * os1], out4);
                  }
                  for (; i0 < (n0+n0_rem); i0++) 
                  {
                    R x0 = I[i0 * is0 + i1 * is1];
                    R x1 = I[i0 * is0 + i1 * is1 + 1];
                    R x2 = I[i0 * is0 + (i1+1) * is1];
                    R x3 = I[i0 * is0 + (i1+1) * is1 + 1];
                    R x4 = I[i0 * is0 + (i1+2) * is1];
                    R x5 = I[i0 * is0 + (i1+2) * is1 + 1];
                    R x6 = I[i0 * is0 + (i1+3) * is1];
                    R x7 = I[i0 * is0 + (i1+3) * is1 + 1];
                    O[i0 * os0 + i1 * os1] = x0;
                    O[i0 * os0 + i1 * os1 + 1] = x1;
                    O[i0 * os0 + (i1+1) * os1] = x2;
                    O[i0 * os0 + (i1+1) * os1 + 1] = x3;
                    O[i0 * os0 + (i1+2) * os1] = x4;
                    O[i0 * os0 + (i1+2) * os1 + 1] = x5;
                    O[i0 * os0 + (i1+3) * os1] = x6;
                    O[i0 * os0 + (i1+3) * os1 + 1] = x7;
                  }
              }
              n0 += n0_rem;
              for (; i1 < (n1+n1_rem); ++i1) 
              {
                  for (i0 = 0; i0 < n0; ++i0) 
                  {
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
              for (i1 = 0; i1 < n1; i1+=4)
              {
                  for (i0 = 0; i0 < n0; i0+=4) 
                  {
                    in1 = _mm256_loadu_ps((float const *)&I[i0 * is0 + i1 * is1]);
                    in2 = _mm256_loadu_ps((float const *)&I[i0 * is0 + (i1+1) * is1]);
                    in3 = _mm256_loadu_ps((float const *)&I[i0 * is0 + (i1+2) * is1]);
                    in4 = _mm256_loadu_ps((float const *)&I[i0 * is0 + (i1+3) * is1]);
    
                    out2 = _mm256_shuffle_ps(in1, in2, 0x44);
                    out4 = _mm256_shuffle_ps(in3, in4, 0x44);
                    out1 = _mm256_permute2f128_ps(out2, out4, 0x20);
                    out3 = _mm256_permute2f128_ps(out2, out4, 0x31);
                    in1 = _mm256_shuffle_ps(in1, in2, 0xEE);
                    in3 = _mm256_shuffle_ps(in3, in4, 0xEE);
                    out2 = _mm256_permute2f128_ps(in1, in3, 0x20);
                    out4 = _mm256_permute2f128_ps(in1, in3, 0x31);
    
                    _mm256_storeu_ps((float *)&O[i0 * os0 + i1 * os1], out1);
                    _mm256_storeu_ps((float *)&O[(i0+1) * os0 + i1 * os1], out2);
                    _mm256_storeu_ps((float *)&O[(i0+2) * os0 + i1 * os1], out3);
                    _mm256_storeu_ps((float *)&O[(i0+3) * os0 + i1 * os1], out4);
                  }
                  for (; i0 < (n0+n0_rem); i0++) 
                  {
                    R x0 = I[i0 * is0 + i1 * is1];
                    R x1 = I[i0 * is0 + i1 * is1 + 1];
                    R x2 = I[i0 * is0 + (i1+1) * is1];
                    R x3 = I[i0 * is0 + (i1+1) * is1 + 1];
                    R x4 = I[i0 * is0 + (i1+2) * is1];
                    R x5 = I[i0 * is0 + (i1+2) * is1 + 1];
                    R x6 = I[i0 * is0 + (i1+3) * is1];
                    R x7 = I[i0 * is0 + (i1+3) * is1 + 1];
                    O[i0 * os0 + i1 * os1] = x0;
                    O[i0 * os0 + i1 * os1 + 1] = x1;
                    O[i0 * os0 + (i1+1) * os1] = x2;
                    O[i0 * os0 + (i1+1) * os1 + 1] = x3;
                    O[i0 * os0 + (i1+2) * os1] = x4;
                    O[i0 * os0 + (i1+2) * os1 + 1] = x5;
                    O[i0 * os0 + (i1+3) * os1] = x6;
                    O[i0 * os0 + (i1+3) * os1 + 1] = x7;
                  }
              }
              n0 += n0_rem;
              for (; i1 < (n1+n1_rem); ++i1) 
              {
                  for (i0 = 0; i0 < n0; ++i0) 
                  {
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
#ifdef AMD_OPT_USE_MEMCPY_TO_CPY
                  memcpy(&O[i1 * os1], &I[i1 * is1], (n0+n0_rem)*2*sizeof(R));
#else
                  for (i0 = 0; i0 < n0; i0+=4) 
                  {
                    in1 = _mm256_loadu_ps((float const *)&I[i0 * is0 + i1 * is1]);
                    _mm256_storeu_ps((float *)&O[i0 * os0 + i1 * os1], in1);
                  }
                  for (; i0 < (n0+n0_rem); i0++) 
                  {
                    R x0 = I[i0 * is0 + i1 * is1];
                    R x1 = I[i0 * is0 + i1 * is1 + 1];
                    O[i0 * os0 + i1 * os1] = x0;
                    O[i0 * os0 + i1 * os1 + 1] = x1;
                  }
#endif
              }
              break;

              default:
              for (i1 = 0; i1 < n1; ++i1)
              {
                for (i0 = 0; i0 < n0; ++i0)
                {
                    R x0 = I0[i0 * is0 + i1 * is1];
                    R x1 = I1[i0 * is0 + i1 * is1];
                    O0[i0 * os0 + i1 * os1] = x0;
                    O1[i0 * os0 + i1 * os1] = x1;
                }
              }
              break;
            }
        }
        break;

        default:
            for (i1 = 0; i1 < n1; ++i1)
            {
                for (i0 = 0; i0 < n0; ++i0)
                {
                    R x0 = I0[i0 * is0 + i1 * is1];
                    R x1 = I1[i0 * is0 + i1 * is1];
                    O0[i0 * os0 + i1 * os1] = x0;
                    O1[i0 * os0 + i1 * os1] = x1;
                }
            }
        break;
     }
}
#else//Default C function
void X(cpy2d_pair)(R *I0, R *I1, R *O0, R *O1,
           INT n0, INT is0, INT os0,
           INT n1, INT is1, INT os1)
{
     INT i0, i1;
     for (i1 = 0; i1 < n1; ++i1)
      for (i0 = 0; i0 < n0; ++i0) {
           R x0 = I0[i0 * is0 + i1 * is1];
           R x1 = I1[i0 * is0 + i1 * is1];
           O0[i0 * os0 + i1 * os1] = x0;
           O1[i0 * os0 + i1 * os1] = x1;
      }
}
#endif//ends

#else//DOUBLE-PRECISION

#if defined(AMD_OPT_IN_PLACE_1D_CPY2D_STABLE_INTRIN)//SIMD optimized function
void X(cpy2d_pair)(R *I0, R *I1, R *O0, R *O1,
           INT n0, INT is0, INT os0,
           INT n1, INT is1, INT os1)
{
     INT i0, i1, v;
     R *I = I0, *O = O0;
     __m256d in1, in2, in3, in4;
     __m256d out1, out2;
     __m128d in1_128, in2_128;
     INT t0, t1, t2, t3;
     INT n0_rem, n1_rem;
     
     switch (((I1 - I0) == 1) & ((O1 - O0) == 1))
     {
        case 1:
        {
            n0_rem = n0&0x1, n1_rem = n1&0x1;
            t0 = (is0==2);
            t1 = (os0==2);
            t2 = (is1==2);
            t3 = (os1==2);
            switch(t0 | (t1 << 1) | (t2 << 2) | (t3 << 3))
            {
            case 1://only is0 is 2. 256-bit contiguous read possible
                n0 = n0 - n0_rem;
                for (i1 = 0; i1 < n1; ++i1) 
                {
                    for (i0 = 0; i0 < n0; i0+=2) 
                    {
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
                for (i1 = 0; i1 < n1; ++i1) 
                {
                    for (i0 = 0; i0 < n0; i0+=2) 
                    {
                        in1_128 = _mm_loadu_pd((double const *)&I[i0 * is0 + i1 * is1]);
                        in2_128 = _mm_loadu_pd((double const *)&I[(i0+1) * is0 + i1 * is1]);
                        in1 = _mm256_castpd128_pd256(in1_128);
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
                    for (i0 = 0; i0 < n0; i0+=2) 
                    {
                        in1 = _mm256_loadu_pd((double const *)&I[i0 * is0 + i1 * is1]);
                        in2 = _mm256_loadu_pd((double const *)&I[(i0+1) * is0 + i1 * is1]);
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
                    for (i0 = 0; i0 < n0; ++i0) 
                    {
                        R x0 = I[i0 * is0 + i1 * is1];
                        R x1 = I[i0 * is0 + i1 * is1 + 1];
                        O[i0 * os0 + i1 * os1] = x0;
                        O[i0 * os0 + i1 * os1 + 1] = x1;
                    }
                }
                break;
    
            case 8://only os1 is 2. 256-bit contiguous write possible
                n1 = n1 - n1_rem;
                for (i1 = 0; i1 < n1; i1+=2) 
                {
                    for (i0 = 0; i0 < n0; ++i0) 
                    {
                        in1_128 = _mm_loadu_pd((double const *)&I[i0 * is0 + i1 * is1]);
                        in2_128 = _mm_loadu_pd((double const *)&I[i0 * is0 + (i1+1) * is1]);
                        in1 = _mm256_castpd128_pd256(in1_128);
                        out1 = _mm256_insertf128_pd(in1, in2_128, 1);
                        _mm256_storeu_pd((double *)&O[i0 * os0 + i1 * os1], out1);
                    }
                }
                if (n1_rem)
                {
                    for (i0 = 0; i0 < n0; ++i0) 
                    {
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
                    for (i0 = 0; i0 < n0; i0+=2) 
                    {
                        in1 = _mm256_loadu_pd((double const *)&I[i0 * is0 + i1 * is1]);
                        in2 = _mm256_loadu_pd((double const *)&I[i0 * is0 + (i1+1) * is1]);
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
                    for (i0 = 0; i0 < n0; ++i0) 
                    {
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
                t1 = n0&0x7;//remainder of 8 in total n0
                t2 = (n0-t1);//multiple of 8
                t3 = (t1)&0x3;//remainder of 4 in remainder of 8
                t1 = t1 - t3;//presence of 4
                t3 = t3 - n0_rem;///presence of 2
                for (i1 = 0; i1 < n1; ++i1)
                {
#ifdef AMD_OPT_USE_MEMCPY_TO_CPY
                    memcpy(&O[i1 * os1], &I[i1 * is1], n0*2*sizeof(R));
#else
                    for (i0 = 0; i0 < t2; i0+=8) 
                    {
                        t0 = i0 * is0 + i1 * is1;
                        in1 = _mm256_loadu_pd((double const *)&I[t0]);
                        t0 += (is0 << 1);
                        in2 = _mm256_loadu_pd((double const *)&I[t0]);
                        t0 += (is0 << 1);
                        in3 = _mm256_loadu_pd((double const *)&I[t0]);
                        t0 += (is0 << 1);
                        in4 = _mm256_loadu_pd((double const *)&I[t0]);
                        t0 = i0 * os0 + i1 * os1;
                        _mm256_storeu_pd((double *)&O[t0], in1);
                        t0 += (os0 << 1);
                        _mm256_storeu_pd((double *)&O[t0], in2);
                        t0 += (os0 << 1);
                        _mm256_storeu_pd((double *)&O[t0], in3);
                        t0 += (os0 << 1);
                        _mm256_storeu_pd((double *)&O[t0], in4);
                    }
                    if (t1) 
                    {
                        t0 = i0 * is0 + i1 * is1;
                        in1 = _mm256_loadu_pd((double const *)&I[t0]);
                        t0 += (is0 << 1);
                        in2 = _mm256_loadu_pd((double const *)&I[t0]);
                        t0 = i0 * os0 + i1 * os1;
                        _mm256_storeu_pd((double *)&O[t0], in1);
                        t0 += (os0 << 1);
                        _mm256_storeu_pd((double *)&O[t0], in2);
                        i0+=4;
                    }
                    if (t3) 
                    {
                        t0 = i0 * is0 + i1 * is1;
                        in1 = _mm256_loadu_pd((double const *)&I[t0]);
                        t0 = i0 * os0 + i1 * os1;
                        _mm256_storeu_pd((double *)&O[t0], in1);
                        i0+=2;
                    }
                    if (n0_rem)
                    {
                        R x0 = I[i0 * is0 + i1 * is1];
                        R x1 = I[i0 * is0 + i1 * is1 + 1];
                        O[i0 * os0 + i1 * os1] = x0;
                        O[i0 * os0 + i1 * os1 + 1] = x1;
                    }
#endif
                }
                break;
    
            default:
                for (i1 = 0; i1 < n1; ++i1)
                {
                    for (i0 = 0; i0 < n0; ++i0)
                    {
                        R x0 = I0[i0 * is0 + i1 * is1];
                        R x1 = I1[i0 * is0 + i1 * is1];
                        O0[i0 * os0 + i1 * os1] = x0;
                        O1[i0 * os0 + i1 * os1] = x1;
                    }
                }
                break;
            }
          }
        break;

        default:
            for (i1 = 0; i1 < n1; ++i1)
            {
                for (i0 = 0; i0 < n0; ++i0)
                {
                    R x0 = I0[i0 * is0 + i1 * is1];
                    R x1 = I1[i0 * is0 + i1 * is1];
                    O0[i0 * os0 + i1 * os1] = x0;
                    O1[i0 * os0 + i1 * os1] = x1;
                }
            }
        break;
     }
}
#else//Default C function
void X(cpy2d_pair)(R *I0, R *I1, R *O0, R *O1,
           INT n0, INT is0, INT os0,
           INT n1, INT is1, INT os1)
{
     INT i0, i1;
     for (i1 = 0; i1 < n1; ++i1)
      for (i0 = 0; i0 < n0; ++i0) {
           R x0 = I0[i0 * is0 + i1 * is1];
           R x1 = I1[i0 * is0 + i1 * is1];
           O0[i0 * os0 + i1 * os1] = x0;
           O1[i0 * os0 + i1 * os1] = x1;
      }
}
#endif//ends
#endif//DOUBLE PRECISION ends

#else //Default(original)

void X(cpy2d_pair)(R *I0, R *I1, R *O0, R *O1,
           INT n0, INT is0, INT os0,
           INT n1, INT is1, INT os1)
{
     INT i0, i1;
     for (i1 = 0; i1 < n1; ++i1)
      for (i0 = 0; i0 < n0; ++i0) {
           R x0 = I0[i0 * is0 + i1 * is1];
           R x1 = I1[i0 * is0 + i1 * is1];
           O0[i0 * os0 + i1 * os1] = x0;
           O1[i0 * os0 + i1 * os1] = x1;
      }
}
#endif

void X(zero1d_pair)(R *O0, R *O1, INT n0, INT os0)
{
     INT i0;
     for (i0 = 0; i0 < n0; ++i0) {
          O0[i0 * os0] = 0;
          O1[i0 * os0] = 0;
     }
}

/* like cpy2d_pair, but read input contiguously if possible */
void X(cpy2d_pair_ci)(R *I0, R *I1, R *O0, R *O1,
              INT n0, INT is0, INT os0,
              INT n1, INT is1, INT os1)
{
     if (IABS(is0) < IABS(is1)) /* inner loop is for n0 */
      X(cpy2d_pair) (I0, I1, O0, O1, n0, is0, os0, n1, is1, os1);
     else
      X(cpy2d_pair) (I0, I1, O0, O1, n1, is1, os1, n0, is0, os0);
}

/* like cpy2d_pair, but write output contiguously if possible */
void X(cpy2d_pair_co)(R *I0, R *I1, R *O0, R *O1,
              INT n0, INT is0, INT os0,
              INT n1, INT is1, INT os1)
{
     if (IABS(os0) < IABS(os1)) /* inner loop is for n0 */
      X(cpy2d_pair) (I0, I1, O0, O1, n0, is0, os0, n1, is1, os1);
     else
      X(cpy2d_pair) (I0, I1, O0, O1, n1, is1, os1, n0, is0, os0);
}
