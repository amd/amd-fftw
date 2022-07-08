/* not worth copyrighting */
/*
 * Copyright (C) 2022, Advanced Micro Devices, Inc. All Rights Reserved.
 *
 */

#include "libbench2/bench.h"

#ifdef AMD_FMV_AUTO
__attribute__((target_clones(TARGET_STRINGS)))
#endif
void caset(bench_complex *A, int n, bench_complex x)
{
     int i;
     for (i = 0; i < n; ++i) {
	  c_re(A[i]) = c_re(x);
	  c_im(A[i]) = c_im(x);
     }
}
