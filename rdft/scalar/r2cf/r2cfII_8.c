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
/* Generated on Thu Dec 10 07:06:02 EST 2020 */

#include "rdft/codelet-rdft.h"

#if defined(ARCH_PREFERS_FMA) || defined(ISA_EXTENSION_PREFERS_FMA)

/* Generated by: ../../../genfft/gen_r2cf.native -fma -compact -variables 4 -pipeline-latency 4 -n 8 -name r2cfII_8 -dft-II -include rdft/scalar/r2cfII.h */

/*
 * This function contains 22 FP additions, 16 FP multiplications,
 * (or, 6 additions, 0 multiplications, 16 fused multiply/add),
 * 18 stack variables, 3 constants, and 16 memory accesses
 */
#include "rdft/scalar/r2cfII.h"

static void r2cfII_8(R *R0, R *R1, R *Cr, R *Ci, stride rs, stride csr, stride csi, INT v, INT ivs, INT ovs)
{
     DK(KP923879532, +0.923879532511286756128183189396788286822416626);
     DK(KP707106781, +0.707106781186547524400844362104849039284835938);
     DK(KP414213562, +0.414213562373095048801688724209698078569671875);
     {
	  INT i;
	  for (i = v; i > 0; i = i - 1, R0 = R0 + ivs, R1 = R1 + ivs, Cr = Cr + ovs, Ci = Ci + ovs, MAKE_VOLATILE_STRIDE(32, rs), MAKE_VOLATILE_STRIDE(32, csr), MAKE_VOLATILE_STRIDE(32, csi)) {
	       E T1, Th, T4, Ti, T8, Te, Tb, Tf, T2, T3;
	       T1 = R0[0];
	       Th = R0[WS(rs, 2)];
	       T2 = R0[WS(rs, 1)];
	       T3 = R0[WS(rs, 3)];
	       T4 = T2 - T3;
	       Ti = T2 + T3;
	       {
		    E T6, T7, T9, Ta;
		    T6 = R1[0];
		    T7 = R1[WS(rs, 2)];
		    T8 = FNMS(KP414213562, T7, T6);
		    Te = FMA(KP414213562, T6, T7);
		    T9 = R1[WS(rs, 3)];
		    Ta = R1[WS(rs, 1)];
		    Tb = FMS(KP414213562, Ta, T9);
		    Tf = FMA(KP414213562, T9, Ta);
	       }
	       {
		    E T5, Tc, Tj, Tk;
		    T5 = FMA(KP707106781, T4, T1);
		    Tc = T8 + Tb;
		    Cr[WS(csr, 3)] = FNMS(KP923879532, Tc, T5);
		    Cr[0] = FMA(KP923879532, Tc, T5);
		    Tj = FMA(KP707106781, Ti, Th);
		    Tk = Te + Tf;
		    Ci[0] = -(FMA(KP923879532, Tk, Tj));
		    Ci[WS(csi, 3)] = FNMS(KP923879532, Tk, Tj);
	       }
	       {
		    E Td, Tg, Tl, Tm;
		    Td = FNMS(KP707106781, T4, T1);
		    Tg = Te - Tf;
		    Cr[WS(csr, 2)] = FNMS(KP923879532, Tg, Td);
		    Cr[WS(csr, 1)] = FMA(KP923879532, Tg, Td);
		    Tl = FNMS(KP707106781, Ti, Th);
		    Tm = Tb - T8;
		    Ci[WS(csi, 2)] = FMS(KP923879532, Tm, Tl);
		    Ci[WS(csi, 1)] = FMA(KP923879532, Tm, Tl);
	       }
	  }
     }
}

static const kr2c_desc desc = { 8, "r2cfII_8", {6, 0, 16, 0}, &GENUS };

void X(codelet_r2cfII_8) (planner *p) { X(kr2c_register) (p, r2cfII_8, &desc);
}

#else

/* Generated by: ../../../genfft/gen_r2cf.native -compact -variables 4 -pipeline-latency 4 -n 8 -name r2cfII_8 -dft-II -include rdft/scalar/r2cfII.h */

/*
 * This function contains 22 FP additions, 10 FP multiplications,
 * (or, 18 additions, 6 multiplications, 4 fused multiply/add),
 * 18 stack variables, 3 constants, and 16 memory accesses
 */
#include "rdft/scalar/r2cfII.h"

static void r2cfII_8(R *R0, R *R1, R *Cr, R *Ci, stride rs, stride csr, stride csi, INT v, INT ivs, INT ovs)
{
     DK(KP382683432, +0.382683432365089771728459984030398866761344562);
     DK(KP923879532, +0.923879532511286756128183189396788286822416626);
     DK(KP707106781, +0.707106781186547524400844362104849039284835938);
     {
	  INT i;
	  for (i = v; i > 0; i = i - 1, R0 = R0 + ivs, R1 = R1 + ivs, Cr = Cr + ovs, Ci = Ci + ovs, MAKE_VOLATILE_STRIDE(32, rs), MAKE_VOLATILE_STRIDE(32, csr), MAKE_VOLATILE_STRIDE(32, csi)) {
	       E T1, Tj, T4, Ti, T8, Te, Tb, Tf, T2, T3;
	       T1 = R0[0];
	       Tj = R0[WS(rs, 2)];
	       T2 = R0[WS(rs, 1)];
	       T3 = R0[WS(rs, 3)];
	       T4 = KP707106781 * (T2 - T3);
	       Ti = KP707106781 * (T2 + T3);
	       {
		    E T6, T7, T9, Ta;
		    T6 = R1[0];
		    T7 = R1[WS(rs, 2)];
		    T8 = FNMS(KP382683432, T7, KP923879532 * T6);
		    Te = FMA(KP382683432, T6, KP923879532 * T7);
		    T9 = R1[WS(rs, 1)];
		    Ta = R1[WS(rs, 3)];
		    Tb = FNMS(KP923879532, Ta, KP382683432 * T9);
		    Tf = FMA(KP923879532, T9, KP382683432 * Ta);
	       }
	       {
		    E T5, Tc, Th, Tk;
		    T5 = T1 + T4;
		    Tc = T8 + Tb;
		    Cr[WS(csr, 3)] = T5 - Tc;
		    Cr[0] = T5 + Tc;
		    Th = Te + Tf;
		    Tk = Ti + Tj;
		    Ci[0] = -(Th + Tk);
		    Ci[WS(csi, 3)] = Tk - Th;
	       }
	       {
		    E Td, Tg, Tl, Tm;
		    Td = T1 - T4;
		    Tg = Te - Tf;
		    Cr[WS(csr, 2)] = Td - Tg;
		    Cr[WS(csr, 1)] = Td + Tg;
		    Tl = Tb - T8;
		    Tm = Tj - Ti;
		    Ci[WS(csi, 2)] = Tl - Tm;
		    Ci[WS(csi, 1)] = Tl + Tm;
	       }
	  }
     }
}

static const kr2c_desc desc = { 8, "r2cfII_8", {18, 6, 4, 0}, &GENUS };

void X(codelet_r2cfII_8) (planner *p) { X(kr2c_register) (p, r2cfII_8, &desc);
}

#endif
