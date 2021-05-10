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
/* Generated on Thu Dec 10 07:06:27 EST 2020 */

#include "rdft/codelet-rdft.h"

#if defined(ARCH_PREFERS_FMA) || defined(ISA_EXTENSION_PREFERS_FMA)

/* Generated by: ../../../genfft/gen_hc2hc.native -fma -compact -variables 4 -pipeline-latency 4 -sign 1 -n 7 -dif -name hb_7 -include rdft/scalar/hb.h */

/*
 * This function contains 72 FP additions, 66 FP multiplications,
 * (or, 18 additions, 12 multiplications, 54 fused multiply/add),
 * 41 stack variables, 6 constants, and 28 memory accesses
 */
#include "rdft/scalar/hb.h"

static void hb_7(R *cr, R *ci, const R *W, stride rs, INT mb, INT me, INT ms)
{
     DK(KP974927912, +0.974927912181823607018131682993931217232785801);
     DK(KP900968867, +0.900968867902419126236102319507445051165919162);
     DK(KP801937735, +0.801937735804838252472204639014890102331838324);
     DK(KP692021471, +0.692021471630095869627814897002069140197260599);
     DK(KP356895867, +0.356895867892209443894399510021300583399127187);
     DK(KP554958132, +0.554958132087371191422194871006410481067288862);
     {
	  INT m;
	  for (m = mb, W = W + ((mb - 1) * 12); m < me; m = m + 1, cr = cr + ms, ci = ci - ms, W = W + 12, MAKE_VOLATILE_STRIDE(14, rs)) {
	       E T1, T4, TC, T7, TB, Ta, TA, TD, TZ, T1l, T1b, TP, Td, Tt, Tw;
	       E Tv, Tu, Tp, Ty, T1j, T1e, TX, TS;
	       T1 = cr[0];
	       {
		    E T2, T3, T1a, TO, Tc;
		    T2 = cr[WS(rs, 1)];
		    T3 = ci[0];
		    T4 = T2 + T3;
		    TC = T2 - T3;
		    {
			 E T5, T6, T8, T9;
			 T5 = cr[WS(rs, 2)];
			 T6 = ci[WS(rs, 1)];
			 T7 = T5 + T6;
			 TB = T5 - T6;
			 T8 = cr[WS(rs, 3)];
			 T9 = ci[WS(rs, 2)];
			 Ta = T8 + T9;
			 TA = T8 - T9;
		    }
		    TD = FNMS(KP554958132, TC, TB);
		    TZ = FMA(KP554958132, TB, TA);
		    T1l = FMA(KP554958132, TA, TC);
		    T1a = FNMS(KP356895867, T7, T4);
		    T1b = FNMS(KP692021471, T1a, Ta);
		    TO = FNMS(KP356895867, T4, Ta);
		    TP = FNMS(KP692021471, TO, T7);
		    Tc = FNMS(KP356895867, Ta, T7);
		    Td = FNMS(KP692021471, Tc, T4);
	       }
	       Tt = ci[WS(rs, 6)];
	       {
		    E Th, Tk, Tn, Tf, Tg;
		    Tf = ci[WS(rs, 3)];
		    Tg = cr[WS(rs, 4)];
		    Th = Tf + Tg;
		    Tw = Tf - Tg;
		    {
			 E Ti, Tj, Tl, Tm;
			 Ti = ci[WS(rs, 4)];
			 Tj = cr[WS(rs, 5)];
			 Tk = Ti + Tj;
			 Tv = Ti - Tj;
			 Tl = ci[WS(rs, 5)];
			 Tm = cr[WS(rs, 6)];
			 Tn = Tl + Tm;
			 Tu = Tl - Tm;
		    }
		    {
			 E To, Tx, T1i, T1d, TW, TR;
			 To = FNMS(KP554958132, Tn, Tk);
			 Tp = FNMS(KP801937735, To, Th);
			 Tx = FNMS(KP356895867, Tw, Tv);
			 Ty = FNMS(KP692021471, Tx, Tu);
			 T1i = FNMS(KP356895867, Tv, Tu);
			 T1j = FNMS(KP692021471, T1i, Tw);
			 T1d = FMA(KP554958132, Th, Tn);
			 T1e = FMA(KP801937735, T1d, Tk);
			 TW = FNMS(KP356895867, Tu, Tw);
			 TX = FNMS(KP692021471, TW, Tv);
			 TR = FMA(KP554958132, Tk, Th);
			 TS = FNMS(KP801937735, TR, Tn);
		    }
	       }
	       cr[0] = T1 + T4 + T7 + Ta;
	       ci[0] = Tt + Tu + Tv + Tw;
	       {
		    E Tq, TI, TF, TL, Te, Tz, TE;
		    Te = FNMS(KP900968867, Td, T1);
		    Tq = FNMS(KP974927912, Tp, Te);
		    TI = FMA(KP974927912, Tp, Te);
		    Tz = FNMS(KP900968867, Ty, Tt);
		    TE = FNMS(KP801937735, TD, TA);
		    TF = FMA(KP974927912, TE, Tz);
		    TL = FNMS(KP974927912, TE, Tz);
		    {
			 E Tb, Tr, Ts, TG;
			 Tb = W[4];
			 Tr = Tb * Tq;
			 Ts = W[5];
			 TG = Ts * Tq;
			 cr[WS(rs, 3)] = FNMS(Ts, TF, Tr);
			 ci[WS(rs, 3)] = FMA(Tb, TF, TG);
		    }
		    {
			 E TH, TJ, TK, TM;
			 TH = W[6];
			 TJ = TH * TI;
			 TK = W[7];
			 TM = TK * TI;
			 cr[WS(rs, 4)] = FNMS(TK, TL, TJ);
			 ci[WS(rs, 4)] = FMA(TH, TL, TM);
		    }
	       }
	       {
		    E TT, T14, T11, T17, TQ, TY, T10;
		    TQ = FNMS(KP900968867, TP, T1);
		    TT = FNMS(KP974927912, TS, TQ);
		    T14 = FMA(KP974927912, TS, TQ);
		    TY = FNMS(KP900968867, TX, Tt);
		    T10 = FNMS(KP801937735, TZ, TC);
		    T11 = FMA(KP974927912, T10, TY);
		    T17 = FNMS(KP974927912, T10, TY);
		    {
			 E TN, TU, TV, T12;
			 TN = W[2];
			 TU = TN * TT;
			 TV = W[3];
			 T12 = TV * TT;
			 cr[WS(rs, 2)] = FNMS(TV, T11, TU);
			 ci[WS(rs, 2)] = FMA(TN, T11, T12);
		    }
		    {
			 E T13, T15, T16, T18;
			 T13 = W[8];
			 T15 = T13 * T14;
			 T16 = W[9];
			 T18 = T16 * T14;
			 cr[WS(rs, 5)] = FNMS(T16, T17, T15);
			 ci[WS(rs, 5)] = FMA(T13, T17, T18);
		    }
	       }
	       {
		    E T1f, T1q, T1n, T1t, T1c, T1k, T1m;
		    T1c = FNMS(KP900968867, T1b, T1);
		    T1f = FNMS(KP974927912, T1e, T1c);
		    T1q = FMA(KP974927912, T1e, T1c);
		    T1k = FNMS(KP900968867, T1j, Tt);
		    T1m = FMA(KP801937735, T1l, TB);
		    T1n = FMA(KP974927912, T1m, T1k);
		    T1t = FNMS(KP974927912, T1m, T1k);
		    {
			 E T19, T1g, T1h, T1o;
			 T19 = W[0];
			 T1g = T19 * T1f;
			 T1h = W[1];
			 T1o = T1h * T1f;
			 cr[WS(rs, 1)] = FNMS(T1h, T1n, T1g);
			 ci[WS(rs, 1)] = FMA(T19, T1n, T1o);
		    }
		    {
			 E T1p, T1r, T1s, T1u;
			 T1p = W[10];
			 T1r = T1p * T1q;
			 T1s = W[11];
			 T1u = T1s * T1q;
			 cr[WS(rs, 6)] = FNMS(T1s, T1t, T1r);
			 ci[WS(rs, 6)] = FMA(T1p, T1t, T1u);
		    }
	       }
	  }
     }
}

static const tw_instr twinstr[] = {
     {TW_FULL, 1, 7},
     {TW_NEXT, 1, 0}
};

static const hc2hc_desc desc = { 7, "hb_7", twinstr, &GENUS, {18, 12, 54, 0} };

void X(codelet_hb_7) (planner *p) {
     X(khc2hc_register) (p, hb_7, &desc);
}
#else

/* Generated by: ../../../genfft/gen_hc2hc.native -compact -variables 4 -pipeline-latency 4 -sign 1 -n 7 -dif -name hb_7 -include rdft/scalar/hb.h */

/*
 * This function contains 72 FP additions, 60 FP multiplications,
 * (or, 36 additions, 24 multiplications, 36 fused multiply/add),
 * 36 stack variables, 6 constants, and 28 memory accesses
 */
#include "rdft/scalar/hb.h"

static void hb_7(R *cr, R *ci, const R *W, stride rs, INT mb, INT me, INT ms)
{
     DK(KP222520933, +0.222520933956314404288902564496794759466355569);
     DK(KP900968867, +0.900968867902419126236102319507445051165919162);
     DK(KP623489801, +0.623489801858733530525004884004239810632274731);
     DK(KP781831482, +0.781831482468029808708444526674057750232334519);
     DK(KP974927912, +0.974927912181823607018131682993931217232785801);
     DK(KP433883739, +0.433883739117558120475768332848358754609990728);
     {
	  INT m;
	  for (m = mb, W = W + ((mb - 1) * 12); m < me; m = m + 1, cr = cr + ms, ci = ci - ms, W = W + 12, MAKE_VOLATILE_STRIDE(14, rs)) {
	       E T1, T4, T7, Ta, Tx, TI, TV, TQ, TE, Tm, Tb, Te, Th, Tk, Tq;
	       E TF, TR, TU, TJ, Tt;
	       {
		    E Tu, Tw, Tv, T2, T3;
		    T1 = cr[0];
		    T2 = cr[WS(rs, 1)];
		    T3 = ci[0];
		    T4 = T2 + T3;
		    Tu = T2 - T3;
		    {
			 E T5, T6, T8, T9;
			 T5 = cr[WS(rs, 2)];
			 T6 = ci[WS(rs, 1)];
			 T7 = T5 + T6;
			 Tw = T5 - T6;
			 T8 = cr[WS(rs, 3)];
			 T9 = ci[WS(rs, 2)];
			 Ta = T8 + T9;
			 Tv = T8 - T9;
		    }
		    Tx = FMA(KP433883739, Tu, KP974927912 * Tv) - (KP781831482 * Tw);
		    TI = FMA(KP781831482, Tu, KP974927912 * Tw) + (KP433883739 * Tv);
		    TV = FNMS(KP781831482, Tv, KP974927912 * Tu) - (KP433883739 * Tw);
		    TQ = FMA(KP623489801, Ta, T1) + FNMA(KP900968867, T7, KP222520933 * T4);
		    TE = FMA(KP623489801, T4, T1) + FNMA(KP900968867, Ta, KP222520933 * T7);
		    Tm = FMA(KP623489801, T7, T1) + FNMA(KP222520933, Ta, KP900968867 * T4);
	       }
	       {
		    E Tp, Tn, To, Tc, Td;
		    Tb = ci[WS(rs, 6)];
		    Tc = ci[WS(rs, 5)];
		    Td = cr[WS(rs, 6)];
		    Te = Tc - Td;
		    Tp = Tc + Td;
		    {
			 E Tf, Tg, Ti, Tj;
			 Tf = ci[WS(rs, 4)];
			 Tg = cr[WS(rs, 5)];
			 Th = Tf - Tg;
			 Tn = Tf + Tg;
			 Ti = ci[WS(rs, 3)];
			 Tj = cr[WS(rs, 4)];
			 Tk = Ti - Tj;
			 To = Ti + Tj;
		    }
		    Tq = FNMS(KP974927912, To, KP781831482 * Tn) - (KP433883739 * Tp);
		    TF = FMA(KP781831482, Tp, KP974927912 * Tn) + (KP433883739 * To);
		    TR = FMA(KP433883739, Tn, KP781831482 * To) - (KP974927912 * Tp);
		    TU = FMA(KP623489801, Tk, Tb) + FNMA(KP900968867, Th, KP222520933 * Te);
		    TJ = FMA(KP623489801, Te, Tb) + FNMA(KP900968867, Tk, KP222520933 * Th);
		    Tt = FMA(KP623489801, Th, Tb) + FNMA(KP222520933, Tk, KP900968867 * Te);
	       }
	       cr[0] = T1 + T4 + T7 + Ta;
	       ci[0] = Tb + Te + Th + Tk;
	       {
		    E Tr, Ty, Tl, Ts;
		    Tr = Tm - Tq;
		    Ty = Tt - Tx;
		    Tl = W[6];
		    Ts = W[7];
		    cr[WS(rs, 4)] = FNMS(Ts, Ty, Tl * Tr);
		    ci[WS(rs, 4)] = FMA(Tl, Ty, Ts * Tr);
	       }
	       {
		    E TY, T10, TX, TZ;
		    TY = TQ + TR;
		    T10 = TV + TU;
		    TX = W[2];
		    TZ = W[3];
		    cr[WS(rs, 2)] = FNMS(TZ, T10, TX * TY);
		    ci[WS(rs, 2)] = FMA(TX, T10, TZ * TY);
	       }
	       {
		    E TA, TC, Tz, TB;
		    TA = Tm + Tq;
		    TC = Tx + Tt;
		    Tz = W[4];
		    TB = W[5];
		    cr[WS(rs, 3)] = FNMS(TB, TC, Tz * TA);
		    ci[WS(rs, 3)] = FMA(Tz, TC, TB * TA);
	       }
	       {
		    E TM, TO, TL, TN;
		    TM = TE + TF;
		    TO = TJ - TI;
		    TL = W[10];
		    TN = W[11];
		    cr[WS(rs, 6)] = FNMS(TN, TO, TL * TM);
		    ci[WS(rs, 6)] = FMA(TL, TO, TN * TM);
	       }
	       {
		    E TS, TW, TP, TT;
		    TS = TQ - TR;
		    TW = TU - TV;
		    TP = W[8];
		    TT = W[9];
		    cr[WS(rs, 5)] = FNMS(TT, TW, TP * TS);
		    ci[WS(rs, 5)] = FMA(TP, TW, TT * TS);
	       }
	       {
		    E TG, TK, TD, TH;
		    TG = TE - TF;
		    TK = TI + TJ;
		    TD = W[0];
		    TH = W[1];
		    cr[WS(rs, 1)] = FNMS(TH, TK, TD * TG);
		    ci[WS(rs, 1)] = FMA(TD, TK, TH * TG);
	       }
	  }
     }
}

static const tw_instr twinstr[] = {
     {TW_FULL, 1, 7},
     {TW_NEXT, 1, 0}
};

static const hc2hc_desc desc = { 7, "hb_7", twinstr, &GENUS, {36, 24, 36, 0} };

void X(codelet_hb_7) (planner *p) {
     X(khc2hc_register) (p, hb_7, &desc);
}
#endif
