/* not worth copyrighting */
/*
 * Copyright (C) 2022, Advanced Micro Devices, Inc. All Rights Reserved.
 *
 */

#include "libbench2/bench.h"

#ifdef AMD_FMV_AUTO
__attribute__((target_clones(TARGET_STRINGS)))
#endif
void aset(bench_real *A, int n, bench_real x)
{
     int i;
     for (i = 0; i < n; ++i)
	  A[i] = x;
}
