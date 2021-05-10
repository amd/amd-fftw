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

/* This file was automatically generated --- DO NOT EDIT */
/* Generated on Thu Dec 10 07:04:41 EST 2020 */

#include "dft/codelet-dft.h"

#if defined(ARCH_PREFERS_FMA) || defined(ISA_EXTENSION_PREFERS_FMA)

/* Generated by: ../../../genfft/gen_notw_c.native -fma -simd -compact -variables 4 -pipeline-latency 8 -n 12 -name n1fv_12 -include dft/simd/n1f.h */

/*
 * This function contains 48 FP additions, 20 FP multiplications,
 * (or, 30 additions, 2 multiplications, 18 fused multiply/add),
 * 27 stack variables, 2 constants, and 24 memory accesses
 */
#include "dft/simd/n1f.h"

static void n1fv_12(const R *ri, const R *ii, R *ro, R *io, stride is, stride os, INT v, INT ivs, INT ovs)
{
     DVK(KP866025403, +0.866025403784438646763723170752936183471402627);
     DVK(KP500000000, +0.500000000000000000000000000000000000000000000);
     {
	  INT i;
	  const R *xi;
	  R *xo;
	  xi = ri;
	  xo = ro;
	  for (i = v; i > 0; i = i - VL, xi = xi + (VL * ivs), xo = xo + (VL * ovs), MAKE_VOLATILE_STRIDE(24, is), MAKE_VOLATILE_STRIDE(24, os)) {
	       V T5, Ta, TG, TF, TB, Tt, Ti, Tm, TJ, TI, TA, Tp;
	       {
		    V T1, T6, T4, Tr, T9, Ts;
		    T1 = LD(&(xi[0]), ivs, &(xi[0]));
		    T6 = LD(&(xi[WS(is, 6)]), ivs, &(xi[0]));
		    {
			 V T2, T3, T7, T8;
			 T2 = LD(&(xi[WS(is, 4)]), ivs, &(xi[0]));
			 T3 = LD(&(xi[WS(is, 8)]), ivs, &(xi[0]));
			 T4 = VADD(T2, T3);
			 Tr = VSUB(T3, T2);
			 T7 = LD(&(xi[WS(is, 10)]), ivs, &(xi[0]));
			 T8 = LD(&(xi[WS(is, 2)]), ivs, &(xi[0]));
			 T9 = VADD(T7, T8);
			 Ts = VSUB(T8, T7);
		    }
		    T5 = VFNMS(LDK(KP500000000), T4, T1);
		    Ta = VFNMS(LDK(KP500000000), T9, T6);
		    TG = VADD(T6, T9);
		    TF = VADD(T1, T4);
		    TB = VADD(Tr, Ts);
		    Tt = VSUB(Tr, Ts);
	       }
	       {
		    V Tk, Tn, Te, Tl, Th, To;
		    Tk = LD(&(xi[WS(is, 3)]), ivs, &(xi[WS(is, 1)]));
		    Tn = LD(&(xi[WS(is, 9)]), ivs, &(xi[WS(is, 1)]));
		    {
			 V Tc, Td, Tf, Tg;
			 Tc = LD(&(xi[WS(is, 11)]), ivs, &(xi[WS(is, 1)]));
			 Td = LD(&(xi[WS(is, 7)]), ivs, &(xi[WS(is, 1)]));
			 Te = VSUB(Tc, Td);
			 Tl = VADD(Td, Tc);
			 Tf = LD(&(xi[WS(is, 1)]), ivs, &(xi[WS(is, 1)]));
			 Tg = LD(&(xi[WS(is, 5)]), ivs, &(xi[WS(is, 1)]));
			 Th = VSUB(Tf, Tg);
			 To = VADD(Tf, Tg);
		    }
		    Ti = VADD(Te, Th);
		    Tm = VFNMS(LDK(KP500000000), Tl, Tk);
		    TJ = VADD(Tn, To);
		    TI = VADD(Tk, Tl);
		    TA = VSUB(Te, Th);
		    Tp = VFNMS(LDK(KP500000000), To, Tn);
	       }
	       {
		    V TH, TK, TL, TM;
		    TH = VSUB(TF, TG);
		    TK = VSUB(TI, TJ);
		    ST(&(xo[WS(os, 9)]), VFNMSI(TK, TH), ovs, &(xo[WS(os, 1)]));
		    ST(&(xo[WS(os, 3)]), VFMAI(TK, TH), ovs, &(xo[WS(os, 1)]));
		    TL = VADD(TF, TG);
		    TM = VADD(TI, TJ);
		    ST(&(xo[WS(os, 6)]), VSUB(TL, TM), ovs, &(xo[0]));
		    ST(&(xo[0]), VADD(TL, TM), ovs, &(xo[0]));
	       }
	       {
		    V Tj, Tv, Tu, Tw, Tb, Tq;
		    Tb = VSUB(T5, Ta);
		    Tj = VFMA(LDK(KP866025403), Ti, Tb);
		    Tv = VFNMS(LDK(KP866025403), Ti, Tb);
		    Tq = VSUB(Tm, Tp);
		    Tu = VFNMS(LDK(KP866025403), Tt, Tq);
		    Tw = VFMA(LDK(KP866025403), Tt, Tq);
		    ST(&(xo[WS(os, 1)]), VFNMSI(Tu, Tj), ovs, &(xo[WS(os, 1)]));
		    ST(&(xo[WS(os, 7)]), VFMAI(Tw, Tv), ovs, &(xo[WS(os, 1)]));
		    ST(&(xo[WS(os, 11)]), VFMAI(Tu, Tj), ovs, &(xo[WS(os, 1)]));
		    ST(&(xo[WS(os, 5)]), VFNMSI(Tw, Tv), ovs, &(xo[WS(os, 1)]));
	       }
	       {
		    V TC, TE, Tz, TD, Tx, Ty;
		    TC = VMUL(LDK(KP866025403), VSUB(TA, TB));
		    TE = VMUL(LDK(KP866025403), VADD(TB, TA));
		    Tx = VADD(T5, Ta);
		    Ty = VADD(Tm, Tp);
		    Tz = VSUB(Tx, Ty);
		    TD = VADD(Tx, Ty);
		    ST(&(xo[WS(os, 2)]), VFMAI(TC, Tz), ovs, &(xo[0]));
		    ST(&(xo[WS(os, 8)]), VFNMSI(TE, TD), ovs, &(xo[0]));
		    ST(&(xo[WS(os, 10)]), VFNMSI(TC, Tz), ovs, &(xo[0]));
		    ST(&(xo[WS(os, 4)]), VFMAI(TE, TD), ovs, &(xo[0]));
	       }
	  }
     }
     VLEAVE();
}

static const kdft_desc desc = { 12, XSIMD_STRING("n1fv_12"), {30, 2, 18, 0}, &GENUS, 0, 0, 0, 0 };

void XSIMD(codelet_n1fv_12) (planner *p) { X(kdft_register) (p, n1fv_12, &desc);
}

#else

/* Generated by: ../../../genfft/gen_notw_c.native -simd -compact -variables 4 -pipeline-latency 8 -n 12 -name n1fv_12 -include dft/simd/n1f.h */

/*
 * This function contains 48 FP additions, 8 FP multiplications,
 * (or, 44 additions, 4 multiplications, 4 fused multiply/add),
 * 27 stack variables, 2 constants, and 24 memory accesses
 */
#include "dft/simd/n1f.h"

static void n1fv_12(const R *ri, const R *ii, R *ro, R *io, stride is, stride os, INT v, INT ivs, INT ovs)
{
     DVK(KP500000000, +0.500000000000000000000000000000000000000000000);
     DVK(KP866025403, +0.866025403784438646763723170752936183471402627);
     {
	  INT i;
	  const R *xi;
	  R *xo;
	  xi = ri;
	  xo = ro;
	  for (i = v; i > 0; i = i - VL, xi = xi + (VL * ivs), xo = xo + (VL * ovs), MAKE_VOLATILE_STRIDE(24, is), MAKE_VOLATILE_STRIDE(24, os)) {
	       V T5, Ta, TJ, Ty, Tq, Tp, Tg, Tl, TI, TA, Tz, Tu;
	       {
		    V T1, T6, T4, Tw, T9, Tx;
		    T1 = LD(&(xi[0]), ivs, &(xi[0]));
		    T6 = LD(&(xi[WS(is, 6)]), ivs, &(xi[0]));
		    {
			 V T2, T3, T7, T8;
			 T2 = LD(&(xi[WS(is, 4)]), ivs, &(xi[0]));
			 T3 = LD(&(xi[WS(is, 8)]), ivs, &(xi[0]));
			 T4 = VADD(T2, T3);
			 Tw = VSUB(T3, T2);
			 T7 = LD(&(xi[WS(is, 10)]), ivs, &(xi[0]));
			 T8 = LD(&(xi[WS(is, 2)]), ivs, &(xi[0]));
			 T9 = VADD(T7, T8);
			 Tx = VSUB(T8, T7);
		    }
		    T5 = VADD(T1, T4);
		    Ta = VADD(T6, T9);
		    TJ = VADD(Tw, Tx);
		    Ty = VMUL(LDK(KP866025403), VSUB(Tw, Tx));
		    Tq = VFNMS(LDK(KP500000000), T9, T6);
		    Tp = VFNMS(LDK(KP500000000), T4, T1);
	       }
	       {
		    V Tc, Th, Tf, Ts, Tk, Tt;
		    Tc = LD(&(xi[WS(is, 3)]), ivs, &(xi[WS(is, 1)]));
		    Th = LD(&(xi[WS(is, 9)]), ivs, &(xi[WS(is, 1)]));
		    {
			 V Td, Te, Ti, Tj;
			 Td = LD(&(xi[WS(is, 7)]), ivs, &(xi[WS(is, 1)]));
			 Te = LD(&(xi[WS(is, 11)]), ivs, &(xi[WS(is, 1)]));
			 Tf = VADD(Td, Te);
			 Ts = VSUB(Te, Td);
			 Ti = LD(&(xi[WS(is, 1)]), ivs, &(xi[WS(is, 1)]));
			 Tj = LD(&(xi[WS(is, 5)]), ivs, &(xi[WS(is, 1)]));
			 Tk = VADD(Ti, Tj);
			 Tt = VSUB(Tj, Ti);
		    }
		    Tg = VADD(Tc, Tf);
		    Tl = VADD(Th, Tk);
		    TI = VADD(Ts, Tt);
		    TA = VFNMS(LDK(KP500000000), Tk, Th);
		    Tz = VFNMS(LDK(KP500000000), Tf, Tc);
		    Tu = VMUL(LDK(KP866025403), VSUB(Ts, Tt));
	       }
	       {
		    V Tb, Tm, Tn, To;
		    Tb = VSUB(T5, Ta);
		    Tm = VBYI(VSUB(Tg, Tl));
		    ST(&(xo[WS(os, 9)]), VSUB(Tb, Tm), ovs, &(xo[WS(os, 1)]));
		    ST(&(xo[WS(os, 3)]), VADD(Tb, Tm), ovs, &(xo[WS(os, 1)]));
		    Tn = VADD(T5, Ta);
		    To = VADD(Tg, Tl);
		    ST(&(xo[WS(os, 6)]), VSUB(Tn, To), ovs, &(xo[0]));
		    ST(&(xo[0]), VADD(Tn, To), ovs, &(xo[0]));
	       }
	       {
		    V Tv, TE, TC, TD, Tr, TB;
		    Tr = VSUB(Tp, Tq);
		    Tv = VSUB(Tr, Tu);
		    TE = VADD(Tr, Tu);
		    TB = VSUB(Tz, TA);
		    TC = VBYI(VADD(Ty, TB));
		    TD = VBYI(VSUB(Ty, TB));
		    ST(&(xo[WS(os, 5)]), VSUB(Tv, TC), ovs, &(xo[WS(os, 1)]));
		    ST(&(xo[WS(os, 11)]), VSUB(TE, TD), ovs, &(xo[WS(os, 1)]));
		    ST(&(xo[WS(os, 7)]), VADD(TC, Tv), ovs, &(xo[WS(os, 1)]));
		    ST(&(xo[WS(os, 1)]), VADD(TD, TE), ovs, &(xo[WS(os, 1)]));
	       }
	       {
		    V TK, TM, TH, TL, TF, TG;
		    TK = VBYI(VMUL(LDK(KP866025403), VSUB(TI, TJ)));
		    TM = VBYI(VMUL(LDK(KP866025403), VADD(TJ, TI)));
		    TF = VADD(Tp, Tq);
		    TG = VADD(Tz, TA);
		    TH = VSUB(TF, TG);
		    TL = VADD(TF, TG);
		    ST(&(xo[WS(os, 10)]), VSUB(TH, TK), ovs, &(xo[0]));
		    ST(&(xo[WS(os, 4)]), VADD(TL, TM), ovs, &(xo[0]));
		    ST(&(xo[WS(os, 2)]), VADD(TH, TK), ovs, &(xo[0]));
		    ST(&(xo[WS(os, 8)]), VSUB(TL, TM), ovs, &(xo[0]));
	       }
	  }
     }
     VLEAVE();
}

static const kdft_desc desc = { 12, XSIMD_STRING("n1fv_12"), {44, 4, 4, 0}, &GENUS, 0, 0, 0, 0 };

void XSIMD(codelet_n1fv_12) (planner *p) { X(kdft_register) (p, n1fv_12, &desc);
}

#endif
