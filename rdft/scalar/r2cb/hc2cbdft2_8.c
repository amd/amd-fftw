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
/* Generated on Thu Dec 10 07:06:50 EST 2020 */

#include "rdft/codelet-rdft.h"

#if defined(ARCH_PREFERS_FMA) || defined(ISA_EXTENSION_PREFERS_FMA)

/* Generated by: ../../../genfft/gen_hc2cdft.native -fma -compact -variables 4 -pipeline-latency 4 -sign 1 -n 8 -dif -name hc2cbdft2_8 -include rdft/scalar/hc2cb.h */

/*
 * This function contains 82 FP additions, 36 FP multiplications,
 * (or, 60 additions, 14 multiplications, 22 fused multiply/add),
 * 41 stack variables, 1 constants, and 32 memory accesses
 */
#include "rdft/scalar/hc2cb.h"

static void hc2cbdft2_8(R *Rp, R *Ip, R *Rm, R *Im, const R *W, stride rs, INT mb, INT me, INT ms)
{
     DK(KP707106781, +0.707106781186547524400844362104849039284835938);
     {
	  INT m;
	  for (m = mb, W = W + ((mb - 1) * 14); m < me; m = m + 1, Rp = Rp + ms, Ip = Ip + ms, Rm = Rm - ms, Im = Im - ms, W = W + 14, MAKE_VOLATILE_STRIDE(32, rs)) {
	       E Tl, T1p, T1g, TM, T1k, TE, TP, T1f, T7, Te, TU, TH, T1l, Tw, T1q;
	       E T1c, T1y;
	       {
		    E T3, TA, Tk, TN, T6, Th, TD, TO, Ta, Tm, Tp, TK, Td, Tr, Tu;
		    E TL, TF, TG;
		    {
			 E T1, T2, Ti, Tj;
			 T1 = Rp[0];
			 T2 = Rm[WS(rs, 3)];
			 T3 = T1 + T2;
			 TA = T1 - T2;
			 Ti = Ip[0];
			 Tj = Im[WS(rs, 3)];
			 Tk = Ti + Tj;
			 TN = Ti - Tj;
		    }
		    {
			 E T4, T5, TB, TC;
			 T4 = Rp[WS(rs, 2)];
			 T5 = Rm[WS(rs, 1)];
			 T6 = T4 + T5;
			 Th = T4 - T5;
			 TB = Ip[WS(rs, 2)];
			 TC = Im[WS(rs, 1)];
			 TD = TB + TC;
			 TO = TB - TC;
		    }
		    {
			 E T8, T9, Tn, To;
			 T8 = Rp[WS(rs, 1)];
			 T9 = Rm[WS(rs, 2)];
			 Ta = T8 + T9;
			 Tm = T8 - T9;
			 Tn = Ip[WS(rs, 1)];
			 To = Im[WS(rs, 2)];
			 Tp = Tn + To;
			 TK = Tn - To;
		    }
		    {
			 E Tb, Tc, Ts, Tt;
			 Tb = Rm[0];
			 Tc = Rp[WS(rs, 3)];
			 Td = Tb + Tc;
			 Tr = Tb - Tc;
			 Ts = Im[0];
			 Tt = Ip[WS(rs, 3)];
			 Tu = Ts + Tt;
			 TL = Tt - Ts;
		    }
		    Tl = Th + Tk;
		    T1p = TA + TD;
		    T1g = TN - TO;
		    TM = TK + TL;
		    T1k = Tk - Th;
		    TE = TA - TD;
		    TP = TN + TO;
		    T1f = Ta - Td;
		    T7 = T3 + T6;
		    Te = Ta + Td;
		    TU = T7 - Te;
		    TF = Tm - Tp;
		    TG = Tr - Tu;
		    TH = TF + TG;
		    T1l = TF - TG;
		    {
			 E Tq, Tv, T1a, T1b;
			 Tq = Tm + Tp;
			 Tv = Tr + Tu;
			 Tw = Tq - Tv;
			 T1q = Tq + Tv;
			 T1a = T3 - T6;
			 T1b = TL - TK;
			 T1c = T1a + T1b;
			 T1y = T1a - T1b;
		    }
	       }
	       {
		    E Tf, TQ, Tx, TI, Ty, TR, Tg, TJ, TS, Tz;
		    Tf = T7 + Te;
		    TQ = TM + TP;
		    Tx = FMA(KP707106781, Tw, Tl);
		    TI = FMA(KP707106781, TH, TE);
		    Tg = W[0];
		    Ty = Tg * Tx;
		    TR = Tg * TI;
		    Tz = W[1];
		    TJ = FMA(Tz, TI, Ty);
		    TS = FNMS(Tz, Tx, TR);
		    Rp[0] = Tf - TJ;
		    Ip[0] = TQ + TS;
		    Rm[0] = Tf + TJ;
		    Im[0] = TS - TQ;
	       }
	       {
		    E T1B, T1A, T1J, T1x, T1z, T1E, T1H, T1F, T1L, T1D;
		    T1B = T1g - T1f;
		    T1A = W[11];
		    T1J = T1A * T1y;
		    T1x = W[10];
		    T1z = T1x * T1y;
		    T1E = FNMS(KP707106781, T1l, T1k);
		    T1H = FMA(KP707106781, T1q, T1p);
		    T1D = W[12];
		    T1F = T1D * T1E;
		    T1L = T1D * T1H;
		    {
			 E T1C, T1K, T1I, T1M, T1G;
			 T1C = FNMS(T1A, T1B, T1z);
			 T1K = FMA(T1x, T1B, T1J);
			 T1G = W[13];
			 T1I = FMA(T1G, T1H, T1F);
			 T1M = FNMS(T1G, T1E, T1L);
			 Rp[WS(rs, 3)] = T1C - T1I;
			 Ip[WS(rs, 3)] = T1K + T1M;
			 Rm[WS(rs, 3)] = T1C + T1I;
			 Im[WS(rs, 3)] = T1M - T1K;
		    }
	       }
	       {
		    E TX, TW, T15, TT, TV, T10, T13, T11, T17, TZ;
		    TX = TP - TM;
		    TW = W[7];
		    T15 = TW * TU;
		    TT = W[6];
		    TV = TT * TU;
		    T10 = FNMS(KP707106781, Tw, Tl);
		    T13 = FNMS(KP707106781, TH, TE);
		    TZ = W[8];
		    T11 = TZ * T10;
		    T17 = TZ * T13;
		    {
			 E TY, T16, T14, T18, T12;
			 TY = FNMS(TW, TX, TV);
			 T16 = FMA(TT, TX, T15);
			 T12 = W[9];
			 T14 = FMA(T12, T13, T11);
			 T18 = FNMS(T12, T10, T17);
			 Rp[WS(rs, 2)] = TY - T14;
			 Ip[WS(rs, 2)] = T16 + T18;
			 Rm[WS(rs, 2)] = TY + T14;
			 Im[WS(rs, 2)] = T18 - T16;
		    }
	       }
	       {
		    E T1h, T1e, T1t, T19, T1d, T1m, T1r, T1n, T1v, T1j;
		    T1h = T1f + T1g;
		    T1e = W[3];
		    T1t = T1e * T1c;
		    T19 = W[2];
		    T1d = T19 * T1c;
		    T1m = FMA(KP707106781, T1l, T1k);
		    T1r = FNMS(KP707106781, T1q, T1p);
		    T1j = W[4];
		    T1n = T1j * T1m;
		    T1v = T1j * T1r;
		    {
			 E T1i, T1u, T1s, T1w, T1o;
			 T1i = FNMS(T1e, T1h, T1d);
			 T1u = FMA(T19, T1h, T1t);
			 T1o = W[5];
			 T1s = FMA(T1o, T1r, T1n);
			 T1w = FNMS(T1o, T1m, T1v);
			 Rp[WS(rs, 1)] = T1i - T1s;
			 Ip[WS(rs, 1)] = T1u + T1w;
			 Rm[WS(rs, 1)] = T1i + T1s;
			 Im[WS(rs, 1)] = T1w - T1u;
		    }
	       }
	  }
     }
}

static const tw_instr twinstr[] = {
     {TW_FULL, 1, 8},
     {TW_NEXT, 1, 0}
};

static const hc2c_desc desc = { 8, "hc2cbdft2_8", twinstr, &GENUS, {60, 14, 22, 0} };

void X(codelet_hc2cbdft2_8) (planner *p) {
     X(khc2c_register) (p, hc2cbdft2_8, &desc, HC2C_VIA_DFT);
}
#else

/* Generated by: ../../../genfft/gen_hc2cdft.native -compact -variables 4 -pipeline-latency 4 -sign 1 -n 8 -dif -name hc2cbdft2_8 -include rdft/scalar/hc2cb.h */

/*
 * This function contains 82 FP additions, 32 FP multiplications,
 * (or, 68 additions, 18 multiplications, 14 fused multiply/add),
 * 30 stack variables, 1 constants, and 32 memory accesses
 */
#include "rdft/scalar/hc2cb.h"

static void hc2cbdft2_8(R *Rp, R *Ip, R *Rm, R *Im, const R *W, stride rs, INT mb, INT me, INT ms)
{
     DK(KP707106781, +0.707106781186547524400844362104849039284835938);
     {
	  INT m;
	  for (m = mb, W = W + ((mb - 1) * 14); m < me; m = m + 1, Rp = Rp + ms, Ip = Ip + ms, Rm = Rm - ms, Im = Im - ms, W = W + 14, MAKE_VOLATILE_STRIDE(32, rs)) {
	       E T7, T1d, T1h, Tl, TG, T14, T19, TO, Te, TL, T18, T15, TB, T1e, Tw;
	       E T1i;
	       {
		    E T3, TC, Tk, TM, T6, Th, TF, TN;
		    {
			 E T1, T2, Ti, Tj;
			 T1 = Rp[0];
			 T2 = Rm[WS(rs, 3)];
			 T3 = T1 + T2;
			 TC = T1 - T2;
			 Ti = Ip[0];
			 Tj = Im[WS(rs, 3)];
			 Tk = Ti + Tj;
			 TM = Ti - Tj;
		    }
		    {
			 E T4, T5, TD, TE;
			 T4 = Rp[WS(rs, 2)];
			 T5 = Rm[WS(rs, 1)];
			 T6 = T4 + T5;
			 Th = T4 - T5;
			 TD = Ip[WS(rs, 2)];
			 TE = Im[WS(rs, 1)];
			 TF = TD + TE;
			 TN = TD - TE;
		    }
		    T7 = T3 + T6;
		    T1d = Tk - Th;
		    T1h = TC + TF;
		    Tl = Th + Tk;
		    TG = TC - TF;
		    T14 = T3 - T6;
		    T19 = TM - TN;
		    TO = TM + TN;
	       }
	       {
		    E Ta, Tm, Tp, TJ, Td, Tr, Tu, TK;
		    {
			 E T8, T9, Tn, To;
			 T8 = Rp[WS(rs, 1)];
			 T9 = Rm[WS(rs, 2)];
			 Ta = T8 + T9;
			 Tm = T8 - T9;
			 Tn = Ip[WS(rs, 1)];
			 To = Im[WS(rs, 2)];
			 Tp = Tn + To;
			 TJ = Tn - To;
		    }
		    {
			 E Tb, Tc, Ts, Tt;
			 Tb = Rm[0];
			 Tc = Rp[WS(rs, 3)];
			 Td = Tb + Tc;
			 Tr = Tb - Tc;
			 Ts = Im[0];
			 Tt = Ip[WS(rs, 3)];
			 Tu = Ts + Tt;
			 TK = Tt - Ts;
		    }
		    Te = Ta + Td;
		    TL = TJ + TK;
		    T18 = Ta - Td;
		    T15 = TK - TJ;
		    {
			 E Tz, TA, Tq, Tv;
			 Tz = Tm - Tp;
			 TA = Tr - Tu;
			 TB = KP707106781 * (Tz + TA);
			 T1e = KP707106781 * (Tz - TA);
			 Tq = Tm + Tp;
			 Tv = Tr + Tu;
			 Tw = KP707106781 * (Tq - Tv);
			 T1i = KP707106781 * (Tq + Tv);
		    }
	       }
	       {
		    E Tf, TP, TI, TQ;
		    Tf = T7 + Te;
		    TP = TL + TO;
		    {
			 E Tx, TH, Tg, Ty;
			 Tx = Tl + Tw;
			 TH = TB + TG;
			 Tg = W[0];
			 Ty = W[1];
			 TI = FMA(Tg, Tx, Ty * TH);
			 TQ = FNMS(Ty, Tx, Tg * TH);
		    }
		    Rp[0] = Tf - TI;
		    Ip[0] = TP + TQ;
		    Rm[0] = Tf + TI;
		    Im[0] = TQ - TP;
	       }
	       {
		    E T1r, T1x, T1w, T1y;
		    {
			 E T1o, T1q, T1n, T1p;
			 T1o = T14 - T15;
			 T1q = T19 - T18;
			 T1n = W[10];
			 T1p = W[11];
			 T1r = FNMS(T1p, T1q, T1n * T1o);
			 T1x = FMA(T1p, T1o, T1n * T1q);
		    }
		    {
			 E T1t, T1v, T1s, T1u;
			 T1t = T1d - T1e;
			 T1v = T1i + T1h;
			 T1s = W[12];
			 T1u = W[13];
			 T1w = FMA(T1s, T1t, T1u * T1v);
			 T1y = FNMS(T1u, T1t, T1s * T1v);
		    }
		    Rp[WS(rs, 3)] = T1r - T1w;
		    Ip[WS(rs, 3)] = T1x + T1y;
		    Rm[WS(rs, 3)] = T1r + T1w;
		    Im[WS(rs, 3)] = T1y - T1x;
	       }
	       {
		    E TV, T11, T10, T12;
		    {
			 E TS, TU, TR, TT;
			 TS = T7 - Te;
			 TU = TO - TL;
			 TR = W[6];
			 TT = W[7];
			 TV = FNMS(TT, TU, TR * TS);
			 T11 = FMA(TT, TS, TR * TU);
		    }
		    {
			 E TX, TZ, TW, TY;
			 TX = Tl - Tw;
			 TZ = TG - TB;
			 TW = W[8];
			 TY = W[9];
			 T10 = FMA(TW, TX, TY * TZ);
			 T12 = FNMS(TY, TX, TW * TZ);
		    }
		    Rp[WS(rs, 2)] = TV - T10;
		    Ip[WS(rs, 2)] = T11 + T12;
		    Rm[WS(rs, 2)] = TV + T10;
		    Im[WS(rs, 2)] = T12 - T11;
	       }
	       {
		    E T1b, T1l, T1k, T1m;
		    {
			 E T16, T1a, T13, T17;
			 T16 = T14 + T15;
			 T1a = T18 + T19;
			 T13 = W[2];
			 T17 = W[3];
			 T1b = FNMS(T17, T1a, T13 * T16);
			 T1l = FMA(T17, T16, T13 * T1a);
		    }
		    {
			 E T1f, T1j, T1c, T1g;
			 T1f = T1d + T1e;
			 T1j = T1h - T1i;
			 T1c = W[4];
			 T1g = W[5];
			 T1k = FMA(T1c, T1f, T1g * T1j);
			 T1m = FNMS(T1g, T1f, T1c * T1j);
		    }
		    Rp[WS(rs, 1)] = T1b - T1k;
		    Ip[WS(rs, 1)] = T1l + T1m;
		    Rm[WS(rs, 1)] = T1b + T1k;
		    Im[WS(rs, 1)] = T1m - T1l;
	       }
	  }
     }
}

static const tw_instr twinstr[] = {
     {TW_FULL, 1, 8},
     {TW_NEXT, 1, 0}
};

static const hc2c_desc desc = { 8, "hc2cbdft2_8", twinstr, &GENUS, {68, 18, 14, 0} };

void X(codelet_hc2cbdft2_8) (planner *p) {
     X(khc2c_register) (p, hc2cbdft2_8, &desc, HC2C_VIA_DFT);
}
#endif
