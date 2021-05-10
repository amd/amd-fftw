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
/* Generated on Thu Dec 10 07:05:57 EST 2020 */

#include "rdft/codelet-rdft.h"

#if defined(ARCH_PREFERS_FMA) || defined(ISA_EXTENSION_PREFERS_FMA)

/* Generated by: ../../../genfft/gen_hc2hc.native -fma -compact -variables 4 -pipeline-latency 4 -twiddle-log3 -precompute-twiddles -n 8 -dit -name hf2_8 -include rdft/scalar/hf.h */

/*
 * This function contains 74 FP additions, 50 FP multiplications,
 * (or, 44 additions, 20 multiplications, 30 fused multiply/add),
 * 48 stack variables, 1 constants, and 32 memory accesses
 */
#include "rdft/scalar/hf.h"

static void hf2_8(R *cr, R *ci, const R *W, stride rs, INT mb, INT me, INT ms)
{
     DK(KP707106781, +0.707106781186547524400844362104849039284835938);
     {
	  INT m;
	  for (m = mb, W = W + ((mb - 1) * 6); m < me; m = m + 1, cr = cr + ms, ci = ci - ms, W = W + 6, MAKE_VOLATILE_STRIDE(16, rs)) {
	       E T2, T3, Tl, Tn, T5, T6, Tf, T7, Ts, Tb, To, Ti, TC, TG;
	       {
		    E T4, Tm, Tr, Ta, TB, TF;
		    T2 = W[0];
		    T3 = W[2];
		    T4 = T2 * T3;
		    Tl = W[4];
		    Tm = T2 * Tl;
		    Tn = W[5];
		    Tr = T2 * Tn;
		    T5 = W[1];
		    T6 = W[3];
		    Ta = T2 * T6;
		    Tf = FMA(T5, T6, T4);
		    T7 = FNMS(T5, T6, T4);
		    Ts = FNMS(T5, Tl, Tr);
		    Tb = FMA(T5, T3, Ta);
		    To = FMA(T5, Tn, Tm);
		    TB = Tf * Tl;
		    TF = Tf * Tn;
		    Ti = FNMS(T5, T3, Ta);
		    TC = FMA(Ti, Tn, TB);
		    TG = FNMS(Ti, Tl, TF);
	       }
	       {
		    E T1, T1s, Td, T1r, Tu, TY, Tk, TW, TN, TR, T18, T1a, T1c, T1d, TA;
		    E TI, T11, T13, T15, T16;
		    T1 = cr[0];
		    T1s = ci[0];
		    {
			 E T8, T9, Tc, T1q;
			 T8 = cr[WS(rs, 4)];
			 T9 = T7 * T8;
			 Tc = ci[WS(rs, 4)];
			 T1q = T7 * Tc;
			 Td = FMA(Tb, Tc, T9);
			 T1r = FNMS(Tb, T8, T1q);
		    }
		    {
			 E Tp, Tq, Tt, TX;
			 Tp = cr[WS(rs, 6)];
			 Tq = To * Tp;
			 Tt = ci[WS(rs, 6)];
			 TX = To * Tt;
			 Tu = FMA(Ts, Tt, Tq);
			 TY = FNMS(Ts, Tp, TX);
		    }
		    {
			 E Tg, Th, Tj, TV;
			 Tg = cr[WS(rs, 2)];
			 Th = Tf * Tg;
			 Tj = ci[WS(rs, 2)];
			 TV = Tf * Tj;
			 Tk = FMA(Ti, Tj, Th);
			 TW = FNMS(Ti, Tg, TV);
		    }
		    {
			 E TK, TL, TM, T19, TO, TP, TQ, T1b;
			 TK = cr[WS(rs, 7)];
			 TL = Tl * TK;
			 TM = ci[WS(rs, 7)];
			 T19 = Tl * TM;
			 TO = cr[WS(rs, 3)];
			 TP = T3 * TO;
			 TQ = ci[WS(rs, 3)];
			 T1b = T3 * TQ;
			 TN = FMA(Tn, TM, TL);
			 TR = FMA(T6, TQ, TP);
			 T18 = TN - TR;
			 T1a = FNMS(Tn, TK, T19);
			 T1c = FNMS(T6, TO, T1b);
			 T1d = T1a - T1c;
		    }
		    {
			 E Tx, Ty, Tz, T12, TD, TE, TH, T14;
			 Tx = cr[WS(rs, 1)];
			 Ty = T2 * Tx;
			 Tz = ci[WS(rs, 1)];
			 T12 = T2 * Tz;
			 TD = cr[WS(rs, 5)];
			 TE = TC * TD;
			 TH = ci[WS(rs, 5)];
			 T14 = TC * TH;
			 TA = FMA(T5, Tz, Ty);
			 TI = FMA(TG, TH, TE);
			 T11 = TA - TI;
			 T13 = FNMS(T5, Tx, T12);
			 T15 = FNMS(TG, TD, T14);
			 T16 = T13 - T15;
		    }
		    {
			 E T10, T1g, T1z, T1B, T1f, T1A, T1j, T1C;
			 {
			      E TU, TZ, T1x, T1y;
			      TU = T1 - Td;
			      TZ = TW - TY;
			      T10 = TU + TZ;
			      T1g = TU - TZ;
			      T1x = Tk - Tu;
			      T1y = T1s - T1r;
			      T1z = T1x + T1y;
			      T1B = T1y - T1x;
			 }
			 {
			      E T17, T1e, T1h, T1i;
			      T17 = T11 + T16;
			      T1e = T18 - T1d;
			      T1f = T17 + T1e;
			      T1A = T1e - T17;
			      T1h = T11 - T16;
			      T1i = T18 + T1d;
			      T1j = T1h + T1i;
			      T1C = T1i - T1h;
			 }
			 ci[WS(rs, 2)] = FNMS(KP707106781, T1f, T10);
			 cr[WS(rs, 5)] = FMS(KP707106781, T1C, T1B);
			 ci[WS(rs, 6)] = FMA(KP707106781, T1C, T1B);
			 cr[WS(rs, 1)] = FMA(KP707106781, T1f, T10);
			 cr[WS(rs, 3)] = FNMS(KP707106781, T1j, T1g);
			 cr[WS(rs, 7)] = FMS(KP707106781, T1A, T1z);
			 ci[WS(rs, 4)] = FMA(KP707106781, T1A, T1z);
			 ci[0] = FMA(KP707106781, T1j, T1g);
		    }
		    {
			 E Tw, T1k, T1u, T1w, TT, T1v, T1n, T1o;
			 {
			      E Te, Tv, T1p, T1t;
			      Te = T1 + Td;
			      Tv = Tk + Tu;
			      Tw = Te + Tv;
			      T1k = Te - Tv;
			      T1p = TW + TY;
			      T1t = T1r + T1s;
			      T1u = T1p + T1t;
			      T1w = T1t - T1p;
			 }
			 {
			      E TJ, TS, T1l, T1m;
			      TJ = TA + TI;
			      TS = TN + TR;
			      TT = TJ + TS;
			      T1v = TS - TJ;
			      T1l = T1a + T1c;
			      T1m = T13 + T15;
			      T1n = T1l - T1m;
			      T1o = T1m + T1l;
			 }
			 ci[WS(rs, 3)] = Tw - TT;
			 cr[WS(rs, 6)] = T1v - T1w;
			 ci[WS(rs, 5)] = T1v + T1w;
			 cr[0] = Tw + TT;
			 cr[WS(rs, 2)] = T1k - T1n;
			 cr[WS(rs, 4)] = T1o - T1u;
			 ci[WS(rs, 7)] = T1o + T1u;
			 ci[WS(rs, 1)] = T1k + T1n;
		    }
	       }
	  }
     }
}

static const tw_instr twinstr[] = {
     {TW_CEXP, 1, 1},
     {TW_CEXP, 1, 3},
     {TW_CEXP, 1, 7},
     {TW_NEXT, 1, 0}
};

static const hc2hc_desc desc = { 8, "hf2_8", twinstr, &GENUS, {44, 20, 30, 0} };

void X(codelet_hf2_8) (planner *p) {
     X(khc2hc_register) (p, hf2_8, &desc);
}
#else

/* Generated by: ../../../genfft/gen_hc2hc.native -compact -variables 4 -pipeline-latency 4 -twiddle-log3 -precompute-twiddles -n 8 -dit -name hf2_8 -include rdft/scalar/hf.h */

/*
 * This function contains 74 FP additions, 44 FP multiplications,
 * (or, 56 additions, 26 multiplications, 18 fused multiply/add),
 * 42 stack variables, 1 constants, and 32 memory accesses
 */
#include "rdft/scalar/hf.h"

static void hf2_8(R *cr, R *ci, const R *W, stride rs, INT mb, INT me, INT ms)
{
     DK(KP707106781, +0.707106781186547524400844362104849039284835938);
     {
	  INT m;
	  for (m = mb, W = W + ((mb - 1) * 6); m < me; m = m + 1, cr = cr + ms, ci = ci - ms, W = W + 6, MAKE_VOLATILE_STRIDE(16, rs)) {
	       E T2, T5, T3, T6, T8, Tc, Tg, Ti, Tl, Tm, Tn, Tz, Tp, Tx;
	       {
		    E T4, Tb, T7, Ta;
		    T2 = W[0];
		    T5 = W[1];
		    T3 = W[2];
		    T6 = W[3];
		    T4 = T2 * T3;
		    Tb = T5 * T3;
		    T7 = T5 * T6;
		    Ta = T2 * T6;
		    T8 = T4 - T7;
		    Tc = Ta + Tb;
		    Tg = T4 + T7;
		    Ti = Ta - Tb;
		    Tl = W[4];
		    Tm = W[5];
		    Tn = FMA(T2, Tl, T5 * Tm);
		    Tz = FNMS(Ti, Tl, Tg * Tm);
		    Tp = FNMS(T5, Tl, T2 * Tm);
		    Tx = FMA(Tg, Tl, Ti * Tm);
	       }
	       {
		    E Tf, T1j, TL, T1d, TJ, T16, TV, TY, Ts, T1i, TO, T1a, TC, T17, TQ;
		    E TT;
		    {
			 E T1, T1c, Te, T1b, T9, Td;
			 T1 = cr[0];
			 T1c = ci[0];
			 T9 = cr[WS(rs, 4)];
			 Td = ci[WS(rs, 4)];
			 Te = FMA(T8, T9, Tc * Td);
			 T1b = FNMS(Tc, T9, T8 * Td);
			 Tf = T1 + Te;
			 T1j = T1c - T1b;
			 TL = T1 - Te;
			 T1d = T1b + T1c;
		    }
		    {
			 E TF, TW, TI, TX;
			 {
			      E TD, TE, TG, TH;
			      TD = cr[WS(rs, 7)];
			      TE = ci[WS(rs, 7)];
			      TF = FMA(Tl, TD, Tm * TE);
			      TW = FNMS(Tm, TD, Tl * TE);
			      TG = cr[WS(rs, 3)];
			      TH = ci[WS(rs, 3)];
			      TI = FMA(T3, TG, T6 * TH);
			      TX = FNMS(T6, TG, T3 * TH);
			 }
			 TJ = TF + TI;
			 T16 = TW + TX;
			 TV = TF - TI;
			 TY = TW - TX;
		    }
		    {
			 E Tk, TM, Tr, TN;
			 {
			      E Th, Tj, To, Tq;
			      Th = cr[WS(rs, 2)];
			      Tj = ci[WS(rs, 2)];
			      Tk = FMA(Tg, Th, Ti * Tj);
			      TM = FNMS(Ti, Th, Tg * Tj);
			      To = cr[WS(rs, 6)];
			      Tq = ci[WS(rs, 6)];
			      Tr = FMA(Tn, To, Tp * Tq);
			      TN = FNMS(Tp, To, Tn * Tq);
			 }
			 Ts = Tk + Tr;
			 T1i = Tk - Tr;
			 TO = TM - TN;
			 T1a = TM + TN;
		    }
		    {
			 E Tw, TR, TB, TS;
			 {
			      E Tu, Tv, Ty, TA;
			      Tu = cr[WS(rs, 1)];
			      Tv = ci[WS(rs, 1)];
			      Tw = FMA(T2, Tu, T5 * Tv);
			      TR = FNMS(T5, Tu, T2 * Tv);
			      Ty = cr[WS(rs, 5)];
			      TA = ci[WS(rs, 5)];
			      TB = FMA(Tx, Ty, Tz * TA);
			      TS = FNMS(Tz, Ty, Tx * TA);
			 }
			 TC = Tw + TB;
			 T17 = TR + TS;
			 TQ = Tw - TB;
			 TT = TR - TS;
		    }
		    {
			 E Tt, TK, T1f, T1g;
			 Tt = Tf + Ts;
			 TK = TC + TJ;
			 ci[WS(rs, 3)] = Tt - TK;
			 cr[0] = Tt + TK;
			 T1f = TJ - TC;
			 T1g = T1d - T1a;
			 cr[WS(rs, 6)] = T1f - T1g;
			 ci[WS(rs, 5)] = T1f + T1g;
			 {
			      E T11, T1m, T14, T1l, T12, T13;
			      T11 = TL - TO;
			      T1m = T1j - T1i;
			      T12 = TQ - TT;
			      T13 = TV + TY;
			      T14 = KP707106781 * (T12 + T13);
			      T1l = KP707106781 * (T13 - T12);
			      cr[WS(rs, 3)] = T11 - T14;
			      ci[WS(rs, 6)] = T1l + T1m;
			      ci[0] = T11 + T14;
			      cr[WS(rs, 5)] = T1l - T1m;
			 }
		    }
		    {
			 E T19, T1e, T15, T18;
			 T19 = T17 + T16;
			 T1e = T1a + T1d;
			 cr[WS(rs, 4)] = T19 - T1e;
			 ci[WS(rs, 7)] = T19 + T1e;
			 T15 = Tf - Ts;
			 T18 = T16 - T17;
			 cr[WS(rs, 2)] = T15 - T18;
			 ci[WS(rs, 1)] = T15 + T18;
			 {
			      E TP, T1k, T10, T1h, TU, TZ;
			      TP = TL + TO;
			      T1k = T1i + T1j;
			      TU = TQ + TT;
			      TZ = TV - TY;
			      T10 = KP707106781 * (TU + TZ);
			      T1h = KP707106781 * (TZ - TU);
			      ci[WS(rs, 2)] = TP - T10;
			      ci[WS(rs, 4)] = T1h + T1k;
			      cr[WS(rs, 1)] = TP + T10;
			      cr[WS(rs, 7)] = T1h - T1k;
			 }
		    }
	       }
	  }
     }
}

static const tw_instr twinstr[] = {
     {TW_CEXP, 1, 1},
     {TW_CEXP, 1, 3},
     {TW_CEXP, 1, 7},
     {TW_NEXT, 1, 0}
};

static const hc2hc_desc desc = { 8, "hf2_8", twinstr, &GENUS, {56, 26, 18, 0} };

void X(codelet_hf2_8) (planner *p) {
     X(khc2hc_register) (p, hf2_8, &desc);
}
#endif
