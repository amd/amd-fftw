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
/* Generated on Thu Dec 10 07:04:11 EST 2020 */

#include "dft/codelet-dft.h"

#if defined(ARCH_PREFERS_FMA) || defined(ISA_EXTENSION_PREFERS_FMA)

/* Generated by: ../../../genfft/gen_twiddle.native -fma -compact -variables 4 -pipeline-latency 4 -n 9 -name t1_9 -include dft/scalar/t.h */

/*
 * This function contains 96 FP additions, 88 FP multiplications,
 * (or, 24 additions, 16 multiplications, 72 fused multiply/add),
 * 55 stack variables, 10 constants, and 36 memory accesses
 */
#include "dft/scalar/t.h"

static void t1_9(R *ri, R *ii, const R *W, stride rs, INT mb, INT me, INT ms)
{
     DK(KP852868531, +0.852868531952443209628250963940074071936020296);
     DK(KP492403876, +0.492403876506104029683371512294761506835321626);
     DK(KP984807753, +0.984807753012208059366743024589523013670643252);
     DK(KP954188894, +0.954188894138671133499268364187245676532219158);
     DK(KP363970234, +0.363970234266202361351047882776834043890471784);
     DK(KP777861913, +0.777861913430206160028177977318626690410586096);
     DK(KP839099631, +0.839099631177280011763127298123181364687434283);
     DK(KP176326980, +0.176326980708464973471090386868618986121633062);
     DK(KP866025403, +0.866025403784438646763723170752936183471402627);
     DK(KP500000000, +0.500000000000000000000000000000000000000000000);
     {
	  INT m;
	  for (m = mb, W = W + (mb * 16); m < me; m = m + 1, ri = ri + ms, ii = ii + ms, W = W + 16, MAKE_VOLATILE_STRIDE(18, rs)) {
	       E T1, T1R, Te, T1W, T10, T1Q, T1l, T1r, Ty, T1p, Tl, T1o, T1g, T1q, T1a;
	       E T1d, TS, T18, TF, T13, T19, T1c;
	       T1 = ri[0];
	       T1R = ii[0];
	       {
		    E T3, T6, T4, TW, T9, Tc, Ta, TY, T2, T8;
		    T3 = ri[WS(rs, 3)];
		    T6 = ii[WS(rs, 3)];
		    T2 = W[4];
		    T4 = T2 * T3;
		    TW = T2 * T6;
		    T9 = ri[WS(rs, 6)];
		    Tc = ii[WS(rs, 6)];
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
			 T1W = Td - T7;
			 T10 = TX - TZ;
			 T1Q = TX + TZ;
		    }
	       }
	       {
		    E Th, Tk, Ti, T1n, Tx, T1i, Tr, T1k, Tg, Tj;
		    Th = ri[WS(rs, 1)];
		    Tk = ii[WS(rs, 1)];
		    Tg = W[0];
		    Ti = Tg * Th;
		    T1n = Tg * Tk;
		    {
			 E Tt, Tw, Tu, T1h, Ts, Tv;
			 Tt = ri[WS(rs, 7)];
			 Tw = ii[WS(rs, 7)];
			 Ts = W[12];
			 Tu = Ts * Tt;
			 T1h = Ts * Tw;
			 Tv = W[13];
			 Tx = FMA(Tv, Tw, Tu);
			 T1i = FNMS(Tv, Tt, T1h);
		    }
		    {
			 E Tn, Tq, To, T1j, Tm, Tp;
			 Tn = ri[WS(rs, 4)];
			 Tq = ii[WS(rs, 4)];
			 Tm = W[6];
			 To = Tm * Tn;
			 T1j = Tm * Tq;
			 Tp = W[7];
			 Tr = FMA(Tp, Tq, To);
			 T1k = FNMS(Tp, Tn, T1j);
		    }
		    T1l = T1i - T1k;
		    T1r = Tr - Tx;
		    Ty = Tr + Tx;
		    T1p = T1k + T1i;
		    Tj = W[1];
		    Tl = FMA(Tj, Tk, Ti);
		    T1o = FNMS(Tj, Th, T1n);
		    T1g = FNMS(KP500000000, Ty, Tl);
		    T1q = FNMS(KP500000000, T1p, T1o);
	       }
	       {
		    E TB, TE, TC, T12, TR, T17, TL, T15, TA, TD;
		    TB = ri[WS(rs, 2)];
		    TE = ii[WS(rs, 2)];
		    TA = W[2];
		    TC = TA * TB;
		    T12 = TA * TE;
		    {
			 E TN, TQ, TO, T16, TM, TP;
			 TN = ri[WS(rs, 8)];
			 TQ = ii[WS(rs, 8)];
			 TM = W[14];
			 TO = TM * TN;
			 T16 = TM * TQ;
			 TP = W[15];
			 TR = FMA(TP, TQ, TO);
			 T17 = FNMS(TP, TN, T16);
		    }
		    {
			 E TH, TK, TI, T14, TG, TJ;
			 TH = ri[WS(rs, 5)];
			 TK = ii[WS(rs, 5)];
			 TG = W[8];
			 TI = TG * TH;
			 T14 = TG * TK;
			 TJ = W[9];
			 TL = FMA(TJ, TK, TI);
			 T15 = FNMS(TJ, TH, T14);
		    }
		    T1a = TR - TL;
		    T1d = T15 - T17;
		    TS = TL + TR;
		    T18 = T15 + T17;
		    TD = W[3];
		    TF = FMA(TD, TE, TC);
		    T13 = FNMS(TD, TB, T12);
		    T19 = FNMS(KP500000000, T18, T13);
		    T1c = FNMS(KP500000000, TS, TF);
	       }
	       {
		    E Tf, T1S, TU, T1U, T1O, T1P, T1L, T1T;
		    Tf = T1 + Te;
		    T1S = T1Q + T1R;
		    {
			 E Tz, TT, T1M, T1N;
			 Tz = Tl + Ty;
			 TT = TF + TS;
			 TU = Tz + TT;
			 T1U = TT - Tz;
			 T1M = T1o + T1p;
			 T1N = T13 + T18;
			 T1O = T1M - T1N;
			 T1P = T1M + T1N;
		    }
		    ri[0] = Tf + TU;
		    ii[0] = T1P + T1S;
		    T1L = FNMS(KP500000000, TU, Tf);
		    ri[WS(rs, 6)] = FNMS(KP866025403, T1O, T1L);
		    ri[WS(rs, 3)] = FMA(KP866025403, T1O, T1L);
		    T1T = FNMS(KP500000000, T1P, T1S);
		    ii[WS(rs, 3)] = FMA(KP866025403, T1U, T1T);
		    ii[WS(rs, 6)] = FNMS(KP866025403, T1U, T1T);
	       }
	       {
		    E T11, T1z, T1X, T21, T1f, T1w, T1t, T1x, T1u, T1Y, T1C, T1I, T1F, T1J, T1G;
		    E T22, TV, T1V;
		    TV = FNMS(KP500000000, Te, T1);
		    T11 = FMA(KP866025403, T10, TV);
		    T1z = FNMS(KP866025403, T10, TV);
		    T1V = FNMS(KP500000000, T1Q, T1R);
		    T1X = FMA(KP866025403, T1W, T1V);
		    T21 = FNMS(KP866025403, T1W, T1V);
		    {
			 E T1b, T1e, T1m, T1s;
			 T1b = FMA(KP866025403, T1a, T19);
			 T1e = FMA(KP866025403, T1d, T1c);
			 T1f = FMA(KP176326980, T1e, T1b);
			 T1w = FNMS(KP176326980, T1b, T1e);
			 T1m = FNMS(KP866025403, T1l, T1g);
			 T1s = FNMS(KP866025403, T1r, T1q);
			 T1t = FMA(KP839099631, T1s, T1m);
			 T1x = FNMS(KP839099631, T1m, T1s);
		    }
		    T1u = FMA(KP777861913, T1t, T1f);
		    T1Y = FNMS(KP777861913, T1x, T1w);
		    {
			 E T1A, T1B, T1D, T1E;
			 T1A = FMA(KP866025403, T1r, T1q);
			 T1B = FMA(KP866025403, T1l, T1g);
			 T1C = FMA(KP176326980, T1B, T1A);
			 T1I = FNMS(KP176326980, T1A, T1B);
			 T1D = FNMS(KP866025403, T1d, T1c);
			 T1E = FNMS(KP866025403, T1a, T19);
			 T1F = FNMS(KP363970234, T1E, T1D);
			 T1J = FMA(KP363970234, T1D, T1E);
		    }
		    T1G = FNMS(KP954188894, T1F, T1C);
		    T22 = FMA(KP954188894, T1J, T1I);
		    ri[WS(rs, 1)] = FMA(KP984807753, T1u, T11);
		    ii[WS(rs, 1)] = FNMS(KP984807753, T1Y, T1X);
		    ri[WS(rs, 2)] = FMA(KP984807753, T1G, T1z);
		    ii[WS(rs, 2)] = FNMS(KP984807753, T22, T21);
		    {
			 E T1v, T1y, T1Z, T20;
			 T1v = FNMS(KP492403876, T1u, T11);
			 T1y = FMA(KP777861913, T1x, T1w);
			 ri[WS(rs, 4)] = FMA(KP852868531, T1y, T1v);
			 ri[WS(rs, 7)] = FNMS(KP852868531, T1y, T1v);
			 T1Z = FMA(KP492403876, T1Y, T1X);
			 T20 = FNMS(KP777861913, T1t, T1f);
			 ii[WS(rs, 4)] = FMA(KP852868531, T20, T1Z);
			 ii[WS(rs, 7)] = FNMS(KP852868531, T20, T1Z);
		    }
		    {
			 E T1H, T1K, T23, T24;
			 T1H = FNMS(KP492403876, T1G, T1z);
			 T1K = FNMS(KP954188894, T1J, T1I);
			 ri[WS(rs, 5)] = FNMS(KP852868531, T1K, T1H);
			 ri[WS(rs, 8)] = FMA(KP852868531, T1K, T1H);
			 T23 = FMA(KP492403876, T22, T21);
			 T24 = FMA(KP954188894, T1F, T1C);
			 ii[WS(rs, 5)] = FNMS(KP852868531, T24, T23);
			 ii[WS(rs, 8)] = FMA(KP852868531, T24, T23);
		    }
	       }
	  }
     }
}

static const tw_instr twinstr[] = {
     {TW_FULL, 0, 9},
     {TW_NEXT, 1, 0}
};

static const ct_desc desc = { 9, "t1_9", twinstr, &GENUS, {24, 16, 72, 0}, 0, 0, 0 };

void X(codelet_t1_9) (planner *p) {
     X(kdft_dit_register) (p, t1_9, &desc);
}
#else

/* Generated by: ../../../genfft/gen_twiddle.native -compact -variables 4 -pipeline-latency 4 -n 9 -name t1_9 -include dft/scalar/t.h */

/*
 * This function contains 96 FP additions, 72 FP multiplications,
 * (or, 60 additions, 36 multiplications, 36 fused multiply/add),
 * 41 stack variables, 8 constants, and 36 memory accesses
 */
#include "dft/scalar/t.h"

static void t1_9(R *ri, R *ii, const R *W, stride rs, INT mb, INT me, INT ms)
{
     DK(KP939692620, +0.939692620785908384054109277324731469936208134);
     DK(KP342020143, +0.342020143325668733044099614682259580763083368);
     DK(KP984807753, +0.984807753012208059366743024589523013670643252);
     DK(KP173648177, +0.173648177666930348851716626769314796000375677);
     DK(KP642787609, +0.642787609686539326322643409907263432907559884);
     DK(KP766044443, +0.766044443118978035202392650555416673935832457);
     DK(KP500000000, +0.500000000000000000000000000000000000000000000);
     DK(KP866025403, +0.866025403784438646763723170752936183471402627);
     {
	  INT m;
	  for (m = mb, W = W + (mb * 16); m < me; m = m + 1, ri = ri + ms, ii = ii + ms, W = W + 16, MAKE_VOLATILE_STRIDE(18, rs)) {
	       E T1, T1B, TQ, T1G, Tc, TN, T1A, T1H, TL, T1x, T17, T1o, T1c, T1n, Tu;
	       E T1w, TW, T1k, T11, T1l;
	       {
		    E T6, TO, Tb, TP;
		    T1 = ri[0];
		    T1B = ii[0];
		    {
			 E T3, T5, T2, T4;
			 T3 = ri[WS(rs, 3)];
			 T5 = ii[WS(rs, 3)];
			 T2 = W[4];
			 T4 = W[5];
			 T6 = FMA(T2, T3, T4 * T5);
			 TO = FNMS(T4, T3, T2 * T5);
		    }
		    {
			 E T8, Ta, T7, T9;
			 T8 = ri[WS(rs, 6)];
			 Ta = ii[WS(rs, 6)];
			 T7 = W[10];
			 T9 = W[11];
			 Tb = FMA(T7, T8, T9 * Ta);
			 TP = FNMS(T9, T8, T7 * Ta);
		    }
		    TQ = KP866025403 * (TO - TP);
		    T1G = KP866025403 * (Tb - T6);
		    Tc = T6 + Tb;
		    TN = FNMS(KP500000000, Tc, T1);
		    T1A = TO + TP;
		    T1H = FNMS(KP500000000, T1A, T1B);
	       }
	       {
		    E Tz, T19, TE, T14, TJ, T15, TK, T1a;
		    {
			 E Tw, Ty, Tv, Tx;
			 Tw = ri[WS(rs, 2)];
			 Ty = ii[WS(rs, 2)];
			 Tv = W[2];
			 Tx = W[3];
			 Tz = FMA(Tv, Tw, Tx * Ty);
			 T19 = FNMS(Tx, Tw, Tv * Ty);
		    }
		    {
			 E TB, TD, TA, TC;
			 TB = ri[WS(rs, 5)];
			 TD = ii[WS(rs, 5)];
			 TA = W[8];
			 TC = W[9];
			 TE = FMA(TA, TB, TC * TD);
			 T14 = FNMS(TC, TB, TA * TD);
		    }
		    {
			 E TG, TI, TF, TH;
			 TG = ri[WS(rs, 8)];
			 TI = ii[WS(rs, 8)];
			 TF = W[14];
			 TH = W[15];
			 TJ = FMA(TF, TG, TH * TI);
			 T15 = FNMS(TH, TG, TF * TI);
		    }
		    TK = TE + TJ;
		    T1a = T14 + T15;
		    TL = Tz + TK;
		    T1x = T19 + T1a;
		    {
			 E T13, T16, T18, T1b;
			 T13 = FNMS(KP500000000, TK, Tz);
			 T16 = KP866025403 * (T14 - T15);
			 T17 = T13 + T16;
			 T1o = T13 - T16;
			 T18 = KP866025403 * (TJ - TE);
			 T1b = FNMS(KP500000000, T1a, T19);
			 T1c = T18 + T1b;
			 T1n = T1b - T18;
		    }
	       }
	       {
		    E Ti, TY, Tn, TT, Ts, TU, Tt, TZ;
		    {
			 E Tf, Th, Te, Tg;
			 Tf = ri[WS(rs, 1)];
			 Th = ii[WS(rs, 1)];
			 Te = W[0];
			 Tg = W[1];
			 Ti = FMA(Te, Tf, Tg * Th);
			 TY = FNMS(Tg, Tf, Te * Th);
		    }
		    {
			 E Tk, Tm, Tj, Tl;
			 Tk = ri[WS(rs, 4)];
			 Tm = ii[WS(rs, 4)];
			 Tj = W[6];
			 Tl = W[7];
			 Tn = FMA(Tj, Tk, Tl * Tm);
			 TT = FNMS(Tl, Tk, Tj * Tm);
		    }
		    {
			 E Tp, Tr, To, Tq;
			 Tp = ri[WS(rs, 7)];
			 Tr = ii[WS(rs, 7)];
			 To = W[12];
			 Tq = W[13];
			 Ts = FMA(To, Tp, Tq * Tr);
			 TU = FNMS(Tq, Tp, To * Tr);
		    }
		    Tt = Tn + Ts;
		    TZ = TT + TU;
		    Tu = Ti + Tt;
		    T1w = TY + TZ;
		    {
			 E TS, TV, TX, T10;
			 TS = FNMS(KP500000000, Tt, Ti);
			 TV = KP866025403 * (TT - TU);
			 TW = TS + TV;
			 T1k = TS - TV;
			 TX = KP866025403 * (Ts - Tn);
			 T10 = FNMS(KP500000000, TZ, TY);
			 T11 = TX + T10;
			 T1l = T10 - TX;
		    }
	       }
	       {
		    E T1y, Td, TM, T1v;
		    T1y = KP866025403 * (T1w - T1x);
		    Td = T1 + Tc;
		    TM = Tu + TL;
		    T1v = FNMS(KP500000000, TM, Td);
		    ri[0] = Td + TM;
		    ri[WS(rs, 3)] = T1v + T1y;
		    ri[WS(rs, 6)] = T1v - T1y;
	       }
	       {
		    E T1D, T1z, T1C, T1E;
		    T1D = KP866025403 * (TL - Tu);
		    T1z = T1w + T1x;
		    T1C = T1A + T1B;
		    T1E = FNMS(KP500000000, T1z, T1C);
		    ii[0] = T1z + T1C;
		    ii[WS(rs, 6)] = T1E - T1D;
		    ii[WS(rs, 3)] = T1D + T1E;
	       }
	       {
		    E TR, T1I, T1e, T1J, T1i, T1F, T1f, T1K;
		    TR = TN + TQ;
		    T1I = T1G + T1H;
		    {
			 E T12, T1d, T1g, T1h;
			 T12 = FMA(KP766044443, TW, KP642787609 * T11);
			 T1d = FMA(KP173648177, T17, KP984807753 * T1c);
			 T1e = T12 + T1d;
			 T1J = KP866025403 * (T1d - T12);
			 T1g = FNMS(KP642787609, TW, KP766044443 * T11);
			 T1h = FNMS(KP984807753, T17, KP173648177 * T1c);
			 T1i = KP866025403 * (T1g - T1h);
			 T1F = T1g + T1h;
		    }
		    ri[WS(rs, 1)] = TR + T1e;
		    ii[WS(rs, 1)] = T1F + T1I;
		    T1f = FNMS(KP500000000, T1e, TR);
		    ri[WS(rs, 7)] = T1f - T1i;
		    ri[WS(rs, 4)] = T1f + T1i;
		    T1K = FNMS(KP500000000, T1F, T1I);
		    ii[WS(rs, 4)] = T1J + T1K;
		    ii[WS(rs, 7)] = T1K - T1J;
	       }
	       {
		    E T1j, T1M, T1q, T1N, T1u, T1L, T1r, T1O;
		    T1j = TN - TQ;
		    T1M = T1H - T1G;
		    {
			 E T1m, T1p, T1s, T1t;
			 T1m = FMA(KP173648177, T1k, KP984807753 * T1l);
			 T1p = FNMS(KP939692620, T1o, KP342020143 * T1n);
			 T1q = T1m + T1p;
			 T1N = KP866025403 * (T1p - T1m);
			 T1s = FNMS(KP984807753, T1k, KP173648177 * T1l);
			 T1t = FMA(KP342020143, T1o, KP939692620 * T1n);
			 T1u = KP866025403 * (T1s + T1t);
			 T1L = T1s - T1t;
		    }
		    ri[WS(rs, 2)] = T1j + T1q;
		    ii[WS(rs, 2)] = T1L + T1M;
		    T1r = FNMS(KP500000000, T1q, T1j);
		    ri[WS(rs, 8)] = T1r - T1u;
		    ri[WS(rs, 5)] = T1r + T1u;
		    T1O = FNMS(KP500000000, T1L, T1M);
		    ii[WS(rs, 5)] = T1N + T1O;
		    ii[WS(rs, 8)] = T1O - T1N;
	       }
	  }
     }
}

static const tw_instr twinstr[] = {
     {TW_FULL, 0, 9},
     {TW_NEXT, 1, 0}
};

static const ct_desc desc = { 9, "t1_9", twinstr, &GENUS, {60, 36, 36, 0}, 0, 0, 0 };

void X(codelet_t1_9) (planner *p) {
     X(kdft_dit_register) (p, t1_9, &desc);
}
#endif
