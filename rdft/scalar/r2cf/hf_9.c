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
/* Generated on Thu Dec 10 07:05:52 EST 2020 */

#include "rdft/codelet-rdft.h"

#if defined(ARCH_PREFERS_FMA) || defined(ISA_EXTENSION_PREFERS_FMA)

/* Generated by: ../../../genfft/gen_hc2hc.native -fma -compact -variables 4 -pipeline-latency 4 -n 9 -dit -name hf_9 -include rdft/scalar/hf.h */

/*
 * This function contains 96 FP additions, 88 FP multiplications,
 * (or, 24 additions, 16 multiplications, 72 fused multiply/add),
 * 55 stack variables, 10 constants, and 36 memory accesses
 */
#include "rdft/scalar/hf.h"

static void hf_9(R *cr, R *ci, const R *W, stride rs, INT mb, INT me, INT ms)
{
     DK(KP852868531, +0.852868531952443209628250963940074071936020296);
     DK(KP492403876, +0.492403876506104029683371512294761506835321626);
     DK(KP984807753, +0.984807753012208059366743024589523013670643252);
     DK(KP777861913, +0.777861913430206160028177977318626690410586096);
     DK(KP839099631, +0.839099631177280011763127298123181364687434283);
     DK(KP954188894, +0.954188894138671133499268364187245676532219158);
     DK(KP363970234, +0.363970234266202361351047882776834043890471784);
     DK(KP176326980, +0.176326980708464973471090386868618986121633062);
     DK(KP866025403, +0.866025403784438646763723170752936183471402627);
     DK(KP500000000, +0.500000000000000000000000000000000000000000000);
     {
	  INT m;
	  for (m = mb, W = W + ((mb - 1) * 16); m < me; m = m + 1, cr = cr + ms, ci = ci - ms, W = W + 16, MAKE_VOLATILE_STRIDE(18, rs)) {
	       E T1, T1P, Te, T1S, T10, T1Q, T1a, T1d, Ty, T18, Tl, T13, T19, T1c, T1l;
	       E T1r, TS, T1p, TF, T1o, T1g, T1q;
	       T1 = cr[0];
	       T1P = ci[0];
	       {
		    E T3, T6, T4, TW, T9, Tc, Ta, TY, T2, T8;
		    T3 = cr[WS(rs, 3)];
		    T6 = ci[WS(rs, 3)];
		    T2 = W[4];
		    T4 = T2 * T3;
		    TW = T2 * T6;
		    T9 = cr[WS(rs, 6)];
		    Tc = ci[WS(rs, 6)];
		    T8 = W[10];
		    Ta = T8 * T9;
		    TY = T8 * Tc;
		    {
			 E T7, TX, Td, TZ, T5, Tb;
			 T5 = W[5];
			 T7 = FMA(T5, T6, T4);
			 TX = FNMS(T5, T3, TW);
			 Tb = W[11];
			 Td = FMA(Tb, Tc, Ta);
			 TZ = FNMS(Tb, T9, TY);
			 Te = T7 + Td;
			 T1S = Td - T7;
			 T10 = TX - TZ;
			 T1Q = TX + TZ;
		    }
	       }
	       {
		    E Th, Tk, Ti, T12, Tx, T17, Tr, T15, Tg, Tj;
		    Th = cr[WS(rs, 1)];
		    Tk = ci[WS(rs, 1)];
		    Tg = W[0];
		    Ti = Tg * Th;
		    T12 = Tg * Tk;
		    {
			 E Tt, Tw, Tu, T16, Ts, Tv;
			 Tt = cr[WS(rs, 7)];
			 Tw = ci[WS(rs, 7)];
			 Ts = W[12];
			 Tu = Ts * Tt;
			 T16 = Ts * Tw;
			 Tv = W[13];
			 Tx = FMA(Tv, Tw, Tu);
			 T17 = FNMS(Tv, Tt, T16);
		    }
		    {
			 E Tn, Tq, To, T14, Tm, Tp;
			 Tn = cr[WS(rs, 4)];
			 Tq = ci[WS(rs, 4)];
			 Tm = W[6];
			 To = Tm * Tn;
			 T14 = Tm * Tq;
			 Tp = W[7];
			 Tr = FMA(Tp, Tq, To);
			 T15 = FNMS(Tp, Tn, T14);
		    }
		    T1a = Tr - Tx;
		    T1d = T15 - T17;
		    Ty = Tr + Tx;
		    T18 = T15 + T17;
		    Tj = W[1];
		    Tl = FMA(Tj, Tk, Ti);
		    T13 = FNMS(Tj, Th, T12);
		    T19 = FNMS(KP500000000, T18, T13);
		    T1c = FNMS(KP500000000, Ty, Tl);
	       }
	       {
		    E TB, TE, TC, T1n, TR, T1k, TL, T1i, TA, TD;
		    TB = cr[WS(rs, 2)];
		    TE = ci[WS(rs, 2)];
		    TA = W[2];
		    TC = TA * TB;
		    T1n = TA * TE;
		    {
			 E TN, TQ, TO, T1j, TM, TP;
			 TN = cr[WS(rs, 8)];
			 TQ = ci[WS(rs, 8)];
			 TM = W[14];
			 TO = TM * TN;
			 T1j = TM * TQ;
			 TP = W[15];
			 TR = FMA(TP, TQ, TO);
			 T1k = FNMS(TP, TN, T1j);
		    }
		    {
			 E TH, TK, TI, T1h, TG, TJ;
			 TH = cr[WS(rs, 5)];
			 TK = ci[WS(rs, 5)];
			 TG = W[8];
			 TI = TG * TH;
			 T1h = TG * TK;
			 TJ = W[9];
			 TL = FMA(TJ, TK, TI);
			 T1i = FNMS(TJ, TH, T1h);
		    }
		    T1l = T1i - T1k;
		    T1r = TR - TL;
		    TS = TL + TR;
		    T1p = T1i + T1k;
		    TD = W[3];
		    TF = FMA(TD, TE, TC);
		    T1o = FNMS(TD, TB, T1n);
		    T1g = FNMS(KP500000000, TS, TF);
		    T1q = FNMS(KP500000000, T1p, T1o);
	       }
	       {
		    E Tf, T21, TU, T24, T1O, T22, T1L, T23;
		    Tf = T1 + Te;
		    T21 = T1Q + T1P;
		    {
			 E Tz, TT, T1M, T1N;
			 Tz = Tl + Ty;
			 TT = TF + TS;
			 TU = Tz + TT;
			 T24 = TT - Tz;
			 T1M = T13 + T18;
			 T1N = T1o + T1p;
			 T1O = T1M - T1N;
			 T22 = T1M + T1N;
		    }
		    cr[0] = Tf + TU;
		    ci[WS(rs, 8)] = T22 + T21;
		    T1L = FNMS(KP500000000, TU, Tf);
		    ci[WS(rs, 2)] = FNMS(KP866025403, T1O, T1L);
		    cr[WS(rs, 3)] = FMA(KP866025403, T1O, T1L);
		    T23 = FNMS(KP500000000, T22, T21);
		    cr[WS(rs, 6)] = FMS(KP866025403, T24, T23);
		    ci[WS(rs, 5)] = FMA(KP866025403, T24, T23);
	       }
	       {
		    E T11, T1z, T1T, T1X, T1f, T1w, T1t, T1x, T1u, T1Y, T1C, T1I, T1F, T1J, T1G;
		    E T1U, TV, T1R;
		    TV = FNMS(KP500000000, Te, T1);
		    T11 = FNMS(KP866025403, T10, TV);
		    T1z = FMA(KP866025403, T10, TV);
		    T1R = FNMS(KP500000000, T1Q, T1P);
		    T1T = FMA(KP866025403, T1S, T1R);
		    T1X = FNMS(KP866025403, T1S, T1R);
		    {
			 E T1b, T1e, T1m, T1s;
			 T1b = FMA(KP866025403, T1a, T19);
			 T1e = FNMS(KP866025403, T1d, T1c);
			 T1f = FMA(KP176326980, T1e, T1b);
			 T1w = FNMS(KP176326980, T1b, T1e);
			 T1m = FNMS(KP866025403, T1l, T1g);
			 T1s = FNMS(KP866025403, T1r, T1q);
			 T1t = FNMS(KP363970234, T1s, T1m);
			 T1x = FMA(KP363970234, T1m, T1s);
		    }
		    T1u = FNMS(KP954188894, T1t, T1f);
		    T1Y = FMA(KP954188894, T1x, T1w);
		    {
			 E T1A, T1B, T1D, T1E;
			 T1A = FMA(KP866025403, T1r, T1q);
			 T1B = FMA(KP866025403, T1l, T1g);
			 T1C = FMA(KP176326980, T1B, T1A);
			 T1I = FNMS(KP176326980, T1A, T1B);
			 T1D = FMA(KP866025403, T1d, T1c);
			 T1E = FNMS(KP866025403, T1a, T19);
			 T1F = FMA(KP839099631, T1E, T1D);
			 T1J = FNMS(KP839099631, T1D, T1E);
		    }
		    T1G = FMA(KP777861913, T1F, T1C);
		    T1U = FNMS(KP777861913, T1J, T1I);
		    cr[WS(rs, 2)] = FMA(KP984807753, T1u, T11);
		    ci[WS(rs, 7)] = FNMS(KP984807753, T1U, T1T);
		    ci[WS(rs, 6)] = FNMS(KP984807753, T1Y, T1X);
		    cr[WS(rs, 1)] = FMA(KP984807753, T1G, T1z);
		    {
			 E T1V, T1W, T1H, T1K;
			 T1V = FMA(KP492403876, T1U, T1T);
			 T1W = FNMS(KP777861913, T1F, T1C);
			 cr[WS(rs, 7)] = FMS(KP852868531, T1W, T1V);
			 ci[WS(rs, 4)] = FMA(KP852868531, T1W, T1V);
			 T1H = FNMS(KP492403876, T1G, T1z);
			 T1K = FMA(KP777861913, T1J, T1I);
			 ci[WS(rs, 1)] = FNMS(KP852868531, T1K, T1H);
			 cr[WS(rs, 4)] = FMA(KP852868531, T1K, T1H);
		    }
		    {
			 E T1v, T1y, T1Z, T20;
			 T1v = FNMS(KP492403876, T1u, T11);
			 T1y = FNMS(KP954188894, T1x, T1w);
			 ci[WS(rs, 3)] = FNMS(KP852868531, T1y, T1v);
			 ci[0] = FMA(KP852868531, T1y, T1v);
			 T1Z = FMA(KP492403876, T1Y, T1X);
			 T20 = FMA(KP954188894, T1t, T1f);
			 cr[WS(rs, 5)] = FMS(KP852868531, T20, T1Z);
			 cr[WS(rs, 8)] = -(FMA(KP852868531, T20, T1Z));
		    }
	       }
	  }
     }
}

static const tw_instr twinstr[] = {
     {TW_FULL, 1, 9},
     {TW_NEXT, 1, 0}
};

static const hc2hc_desc desc = { 9, "hf_9", twinstr, &GENUS, {24, 16, 72, 0} };

void X(codelet_hf_9) (planner *p) {
     X(khc2hc_register) (p, hf_9, &desc);
}
#else

/* Generated by: ../../../genfft/gen_hc2hc.native -compact -variables 4 -pipeline-latency 4 -n 9 -dit -name hf_9 -include rdft/scalar/hf.h */

/*
 * This function contains 96 FP additions, 72 FP multiplications,
 * (or, 60 additions, 36 multiplications, 36 fused multiply/add),
 * 41 stack variables, 8 constants, and 36 memory accesses
 */
#include "rdft/scalar/hf.h"

static void hf_9(R *cr, R *ci, const R *W, stride rs, INT mb, INT me, INT ms)
{
     DK(KP642787609, +0.642787609686539326322643409907263432907559884);
     DK(KP766044443, +0.766044443118978035202392650555416673935832457);
     DK(KP939692620, +0.939692620785908384054109277324731469936208134);
     DK(KP342020143, +0.342020143325668733044099614682259580763083368);
     DK(KP984807753, +0.984807753012208059366743024589523013670643252);
     DK(KP173648177, +0.173648177666930348851716626769314796000375677);
     DK(KP500000000, +0.500000000000000000000000000000000000000000000);
     DK(KP866025403, +0.866025403784438646763723170752936183471402627);
     {
	  INT m;
	  for (m = mb, W = W + ((mb - 1) * 16); m < me; m = m + 1, cr = cr + ms, ci = ci - ms, W = W + 16, MAKE_VOLATILE_STRIDE(18, rs)) {
	       E T1, T1B, TQ, T1A, Tc, TN, T1C, T1D, TL, T1x, T19, T1o, T1c, T1n, Tu;
	       E T1w, TW, T1k, T11, T1l;
	       {
		    E T6, TO, Tb, TP;
		    T1 = cr[0];
		    T1B = ci[0];
		    {
			 E T3, T5, T2, T4;
			 T3 = cr[WS(rs, 3)];
			 T5 = ci[WS(rs, 3)];
			 T2 = W[4];
			 T4 = W[5];
			 T6 = FMA(T2, T3, T4 * T5);
			 TO = FNMS(T4, T3, T2 * T5);
		    }
		    {
			 E T8, Ta, T7, T9;
			 T8 = cr[WS(rs, 6)];
			 Ta = ci[WS(rs, 6)];
			 T7 = W[10];
			 T9 = W[11];
			 Tb = FMA(T7, T8, T9 * Ta);
			 TP = FNMS(T9, T8, T7 * Ta);
		    }
		    TQ = KP866025403 * (TO - TP);
		    T1A = KP866025403 * (Tb - T6);
		    Tc = T6 + Tb;
		    TN = FNMS(KP500000000, Tc, T1);
		    T1C = TO + TP;
		    T1D = FNMS(KP500000000, T1C, T1B);
	       }
	       {
		    E Tz, T13, TE, T14, TJ, T15, TK, T16;
		    {
			 E Tw, Ty, Tv, Tx;
			 Tw = cr[WS(rs, 2)];
			 Ty = ci[WS(rs, 2)];
			 Tv = W[2];
			 Tx = W[3];
			 Tz = FMA(Tv, Tw, Tx * Ty);
			 T13 = FNMS(Tx, Tw, Tv * Ty);
		    }
		    {
			 E TB, TD, TA, TC;
			 TB = cr[WS(rs, 5)];
			 TD = ci[WS(rs, 5)];
			 TA = W[8];
			 TC = W[9];
			 TE = FMA(TA, TB, TC * TD);
			 T14 = FNMS(TC, TB, TA * TD);
		    }
		    {
			 E TG, TI, TF, TH;
			 TG = cr[WS(rs, 8)];
			 TI = ci[WS(rs, 8)];
			 TF = W[14];
			 TH = W[15];
			 TJ = FMA(TF, TG, TH * TI);
			 T15 = FNMS(TH, TG, TF * TI);
		    }
		    TK = TE + TJ;
		    T16 = T14 + T15;
		    TL = Tz + TK;
		    T1x = T13 + T16;
		    {
			 E T17, T18, T1a, T1b;
			 T17 = FNMS(KP500000000, T16, T13);
			 T18 = KP866025403 * (TJ - TE);
			 T19 = T17 - T18;
			 T1o = T18 + T17;
			 T1a = FNMS(KP500000000, TK, Tz);
			 T1b = KP866025403 * (T14 - T15);
			 T1c = T1a - T1b;
			 T1n = T1a + T1b;
		    }
	       }
	       {
		    E Ti, TX, Tn, TT, Ts, TU, Tt, TY;
		    {
			 E Tf, Th, Te, Tg;
			 Tf = cr[WS(rs, 1)];
			 Th = ci[WS(rs, 1)];
			 Te = W[0];
			 Tg = W[1];
			 Ti = FMA(Te, Tf, Tg * Th);
			 TX = FNMS(Tg, Tf, Te * Th);
		    }
		    {
			 E Tk, Tm, Tj, Tl;
			 Tk = cr[WS(rs, 4)];
			 Tm = ci[WS(rs, 4)];
			 Tj = W[6];
			 Tl = W[7];
			 Tn = FMA(Tj, Tk, Tl * Tm);
			 TT = FNMS(Tl, Tk, Tj * Tm);
		    }
		    {
			 E Tp, Tr, To, Tq;
			 Tp = cr[WS(rs, 7)];
			 Tr = ci[WS(rs, 7)];
			 To = W[12];
			 Tq = W[13];
			 Ts = FMA(To, Tp, Tq * Tr);
			 TU = FNMS(Tq, Tp, To * Tr);
		    }
		    Tt = Tn + Ts;
		    TY = TT + TU;
		    Tu = Ti + Tt;
		    T1w = TX + TY;
		    {
			 E TS, TV, TZ, T10;
			 TS = FNMS(KP500000000, Tt, Ti);
			 TV = KP866025403 * (TT - TU);
			 TW = TS - TV;
			 T1k = TS + TV;
			 TZ = FNMS(KP500000000, TY, TX);
			 T10 = KP866025403 * (Ts - Tn);
			 T11 = TZ - T10;
			 T1l = T10 + TZ;
		    }
	       }
	       {
		    E T1y, Td, TM, T1v;
		    T1y = KP866025403 * (T1w - T1x);
		    Td = T1 + Tc;
		    TM = Tu + TL;
		    T1v = FNMS(KP500000000, TM, Td);
		    cr[0] = Td + TM;
		    cr[WS(rs, 3)] = T1v + T1y;
		    ci[WS(rs, 2)] = T1v - T1y;
	       }
	       {
		    E TR, T1I, T1e, T1K, T1i, T1H, T1f, T1J;
		    TR = TN - TQ;
		    T1I = T1D - T1A;
		    {
			 E T12, T1d, T1g, T1h;
			 T12 = FMA(KP173648177, TW, KP984807753 * T11);
			 T1d = FNMS(KP939692620, T1c, KP342020143 * T19);
			 T1e = T12 + T1d;
			 T1K = KP866025403 * (T1d - T12);
			 T1g = FNMS(KP984807753, TW, KP173648177 * T11);
			 T1h = FMA(KP342020143, T1c, KP939692620 * T19);
			 T1i = KP866025403 * (T1g + T1h);
			 T1H = T1g - T1h;
		    }
		    cr[WS(rs, 2)] = TR + T1e;
		    ci[WS(rs, 6)] = T1H + T1I;
		    T1f = FNMS(KP500000000, T1e, TR);
		    ci[0] = T1f - T1i;
		    ci[WS(rs, 3)] = T1f + T1i;
		    T1J = FMS(KP500000000, T1H, T1I);
		    cr[WS(rs, 5)] = T1J - T1K;
		    cr[WS(rs, 8)] = T1K + T1J;
	       }
	       {
		    E T1L, T1M, T1N, T1O;
		    T1L = KP866025403 * (TL - Tu);
		    T1M = T1C + T1B;
		    T1N = T1w + T1x;
		    T1O = FNMS(KP500000000, T1N, T1M);
		    cr[WS(rs, 6)] = T1L - T1O;
		    ci[WS(rs, 8)] = T1N + T1M;
		    ci[WS(rs, 5)] = T1L + T1O;
	       }
	       {
		    E T1j, T1E, T1q, T1z, T1u, T1F, T1r, T1G;
		    T1j = TN + TQ;
		    T1E = T1A + T1D;
		    {
			 E T1m, T1p, T1s, T1t;
			 T1m = FMA(KP766044443, T1k, KP642787609 * T1l);
			 T1p = FMA(KP173648177, T1n, KP984807753 * T1o);
			 T1q = T1m + T1p;
			 T1z = KP866025403 * (T1p - T1m);
			 T1s = FNMS(KP642787609, T1k, KP766044443 * T1l);
			 T1t = FNMS(KP984807753, T1n, KP173648177 * T1o);
			 T1u = KP866025403 * (T1s - T1t);
			 T1F = T1s + T1t;
		    }
		    cr[WS(rs, 1)] = T1j + T1q;
		    T1r = FNMS(KP500000000, T1q, T1j);
		    ci[WS(rs, 1)] = T1r - T1u;
		    cr[WS(rs, 4)] = T1r + T1u;
		    ci[WS(rs, 7)] = T1F + T1E;
		    T1G = FNMS(KP500000000, T1F, T1E);
		    cr[WS(rs, 7)] = T1z - T1G;
		    ci[WS(rs, 4)] = T1z + T1G;
	       }
	  }
     }
}

static const tw_instr twinstr[] = {
     {TW_FULL, 1, 9},
     {TW_NEXT, 1, 0}
};

static const hc2hc_desc desc = { 9, "hf_9", twinstr, &GENUS, {60, 36, 36, 0} };

void X(codelet_hf_9) (planner *p) {
     X(khc2hc_register) (p, hf_9, &desc);
}
#endif
