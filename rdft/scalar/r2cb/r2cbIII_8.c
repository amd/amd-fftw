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
/* Generated on Thu Dec 10 07:06:37 EST 2020 */

#include "rdft/codelet-rdft.h"

#if defined(ARCH_PREFERS_FMA) || defined(ISA_EXTENSION_PREFERS_FMA)

/* Generated by: ../../../genfft/gen_r2cb.native -fma -compact -variables 4 -pipeline-latency 4 -sign 1 -n 8 -name r2cbIII_8 -dft-III -include rdft/scalar/r2cbIII.h */

/*
 * This function contains 22 FP additions, 12 FP multiplications,
 * (or, 18 additions, 8 multiplications, 4 fused multiply/add),
 * 19 stack variables, 4 constants, and 16 memory accesses
 */
#include "rdft/scalar/r2cbIII.h"

static void r2cbIII_8(R *R0, R *R1, R *Cr, R *Ci, stride rs, stride csr, stride csi, INT v, INT ivs, INT ovs)
{
     DK(KP1_414213562, +1.414213562373095048801688724209698078569671875);
     DK(KP414213562, +0.414213562373095048801688724209698078569671875);
     DK(KP1_847759065, +1.847759065022573512256366378793576573644833252);
     DK(KP2_000000000, +2.000000000000000000000000000000000000000000000);
     {
	  INT i;
	  for (i = v; i > 0; i = i - 1, R0 = R0 + ovs, R1 = R1 + ovs, Cr = Cr + ivs, Ci = Ci + ivs, MAKE_VOLATILE_STRIDE(32, rs), MAKE_VOLATILE_STRIDE(32, csr), MAKE_VOLATILE_STRIDE(32, csi)) {
	       E T3, T7, Tf, Tl, T6, Tc, Ta, Tk, Tb, Tg;
	       {
		    E T1, T2, Td, Te;
		    T1 = Cr[0];
		    T2 = Cr[WS(csr, 3)];
		    T3 = T1 + T2;
		    T7 = T1 - T2;
		    Td = Ci[0];
		    Te = Ci[WS(csi, 3)];
		    Tf = Td + Te;
		    Tl = Te - Td;
	       }
	       {
		    E T4, T5, T8, T9;
		    T4 = Cr[WS(csr, 2)];
		    T5 = Cr[WS(csr, 1)];
		    T6 = T4 + T5;
		    Tc = T4 - T5;
		    T8 = Ci[WS(csi, 2)];
		    T9 = Ci[WS(csi, 1)];
		    Ta = T8 + T9;
		    Tk = T8 - T9;
	       }
	       R0[0] = KP2_000000000 * (T3 + T6);
	       R0[WS(rs, 2)] = KP2_000000000 * (Tl - Tk);
	       Tb = T7 - Ta;
	       Tg = Tc + Tf;
	       R1[0] = KP1_847759065 * (FNMS(KP414213562, Tg, Tb));
	       R1[WS(rs, 2)] = -(KP1_847759065 * (FMA(KP414213562, Tb, Tg)));
	       {
		    E Th, Ti, Tj, Tm;
		    Th = Tc - Tf;
		    Ti = T7 + Ta;
		    R1[WS(rs, 1)] = KP1_847759065 * (FMA(KP414213562, Ti, Th));
		    R1[WS(rs, 3)] = -(KP1_847759065 * (FNMS(KP414213562, Th, Ti)));
		    Tj = T3 - T6;
		    Tm = Tk + Tl;
		    R0[WS(rs, 1)] = KP1_414213562 * (Tj + Tm);
		    R0[WS(rs, 3)] = KP1_414213562 * (Tm - Tj);
	       }
	  }
     }
}

static const kr2c_desc desc = { 8, "r2cbIII_8", {18, 8, 4, 0}, &GENUS };

void X(codelet_r2cbIII_8) (planner *p) { X(kr2c_register) (p, r2cbIII_8, &desc);
}

#else

/* Generated by: ../../../genfft/gen_r2cb.native -compact -variables 4 -pipeline-latency 4 -sign 1 -n 8 -name r2cbIII_8 -dft-III -include rdft/scalar/r2cbIII.h */

/*
 * This function contains 22 FP additions, 12 FP multiplications,
 * (or, 18 additions, 8 multiplications, 4 fused multiply/add),
 * 19 stack variables, 4 constants, and 16 memory accesses
 */
#include "rdft/scalar/r2cbIII.h"

static void r2cbIII_8(R *R0, R *R1, R *Cr, R *Ci, stride rs, stride csr, stride csi, INT v, INT ivs, INT ovs)
{
     DK(KP1_414213562, +1.414213562373095048801688724209698078569671875);
     DK(KP765366864, +0.765366864730179543456919968060797733522689125);
     DK(KP1_847759065, +1.847759065022573512256366378793576573644833252);
     DK(KP2_000000000, +2.000000000000000000000000000000000000000000000);
     {
	  INT i;
	  for (i = v; i > 0; i = i - 1, R0 = R0 + ovs, R1 = R1 + ovs, Cr = Cr + ivs, Ci = Ci + ivs, MAKE_VOLATILE_STRIDE(32, rs), MAKE_VOLATILE_STRIDE(32, csr), MAKE_VOLATILE_STRIDE(32, csi)) {
	       E T3, T7, Tf, Tl, T6, Tc, Ta, Tk, Tb, Tg;
	       {
		    E T1, T2, Td, Te;
		    T1 = Cr[0];
		    T2 = Cr[WS(csr, 3)];
		    T3 = T1 + T2;
		    T7 = T1 - T2;
		    Td = Ci[0];
		    Te = Ci[WS(csi, 3)];
		    Tf = Td + Te;
		    Tl = Te - Td;
	       }
	       {
		    E T4, T5, T8, T9;
		    T4 = Cr[WS(csr, 2)];
		    T5 = Cr[WS(csr, 1)];
		    T6 = T4 + T5;
		    Tc = T4 - T5;
		    T8 = Ci[WS(csi, 2)];
		    T9 = Ci[WS(csi, 1)];
		    Ta = T8 + T9;
		    Tk = T8 - T9;
	       }
	       R0[0] = KP2_000000000 * (T3 + T6);
	       R0[WS(rs, 2)] = KP2_000000000 * (Tl - Tk);
	       Tb = T7 - Ta;
	       Tg = Tc + Tf;
	       R1[0] = FNMS(KP765366864, Tg, KP1_847759065 * Tb);
	       R1[WS(rs, 2)] = -(FMA(KP765366864, Tb, KP1_847759065 * Tg));
	       {
		    E Th, Ti, Tj, Tm;
		    Th = T7 + Ta;
		    Ti = Tc - Tf;
		    R1[WS(rs, 1)] = FMA(KP765366864, Th, KP1_847759065 * Ti);
		    R1[WS(rs, 3)] = FNMS(KP1_847759065, Th, KP765366864 * Ti);
		    Tj = T3 - T6;
		    Tm = Tk + Tl;
		    R0[WS(rs, 1)] = KP1_414213562 * (Tj + Tm);
		    R0[WS(rs, 3)] = KP1_414213562 * (Tm - Tj);
	       }
	  }
     }
}

static const kr2c_desc desc = { 8, "r2cbIII_8", {18, 8, 4, 0}, &GENUS };

void X(codelet_r2cbIII_8) (planner *p) { X(kr2c_register) (p, r2cbIII_8, &desc);
}

#endif
