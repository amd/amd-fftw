cmake_minimum_required (VERSION 3.15)

if (NOT DEFINED CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE Release CACHE STRING "Build type")
endif ()

project (fftw)

if (POLICY CMP0042)
  cmake_policy (SET CMP0042 NEW)
endif ()

SET(AMD_ARCH "znver1" CACHE STRING "select AMD zen version for Clang toolchain")

if (CMAKE_C_COMPILER_ID MATCHES Clang)
  if ("${AMD_ARCH}" STREQUAL "")
    message(FATAL_ERROR "Machine arch missing! Select one of znver1, znver2, znver3 or znver4")
  elseif (${AMD_ARCH} STREQUAL "znver1")
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=znver1")
  elseif (${AMD_ARCH} STREQUAL "znver2")
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=znver2")
  elseif (${AMD_ARCH} STREQUAL "znver3")
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=znver3")
  elseif (${AMD_ARCH} STREQUAL "znver4")
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=znver4")
  else ()
    message(FATAL_ERROR "Unsupported Machine arch! Select one of znver1, znver2, znver3 or znver4")
  endif ()
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mno-prfchw")
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mno-prfchw")
endif ()

option (BUILD_SHARED_LIBS "Build shared libraries" ON)
option (BUILD_TESTS "Build tests" ON)
option (ENABLE_VERBOSE_MODE "Enable verbose mode" ON)

option (ENABLE_OPENMP "Use OpenMP for multithreading" OFF)
option (ENABLE_THREADS "Use threads for multithreading" OFF)
option (WITH_COMBINED_THREADS "Merge thread library" OFF)

option (ENABLE_FLOAT "single-precision" OFF)
option (ENABLE_LONG_DOUBLE "long-double precision" OFF)
# Not available in Windows due to limitations on 128 bit datatype
#option (ENABLE_QUAD_PRECISION "quadruple-precision" OFF)

option (ENABLE_SSE "Compile with SSE instruction set support" OFF)
option (ENABLE_SSE2 "Compile with SSE2 instruction set support" OFF)
option (ENABLE_AVX "Compile with AVX instruction set support" OFF)
option (ENABLE_AVX2 "Compile with AVX2 instruction set support" ON)
option (ENABLE_AVX512 "Compile with AVX512 instruction set support" OFF)

option (ENABLE_AMD_OPT "Enable AMD specific optimization" OFF)
option (ENABLE_AMD_MPIFFT "Enable AMD MPI FFT optimizations" OFF)
option (ENABLE_AMD_MPI_VADER_LIMIT "Enable advanced AMD MPI optimizations for in-place transpose based on VADER LIMIT" OFF)
option (ENABLE_MPI "compile FFTW MPI library" OFF)

option (ENABLE_AMD_TRANS "Enable amd optimized transpose" OFF)

option (ENABLE_AMD_FAST_PLANNER "Enable for a faster planning time on AMD cpus" OFF)
option (ENABLE_AMD_TOP_N_planner "Enable AMD Top N Planner for AMD cpus" OFF)

option (ENABLE_AMD_APP_OPT "Enable AMD application optimization layer for HPC and scientific applications" OFF)

if(ENABLE_VERBOSE_MODE)
  if(CMAKE_C_COMPILER_ID MATCHES MSVC)
    set(CMAKE_VERBOSE_MAKEFILE ON CACHE BOOL "ON" FORCE)
  elseif(CMAKE_C_COMPILER_ID MATCHES Clang)
    add_compile_options(-v)
  endif()
endif()

# Release mode setting for optimization flags in MSVC and Clang (in Debug mode optimization is disabled)
if (CMAKE_C_COMPILER_ID MATCHES MSVC OR CMAKE_C_COMPILER_ID MATCHES Clang)
  if(CMAKE_BUILD_TYPE STREQUAL "Release")
    set(CMAKE_CXX_FLAGS_RELEASE "/O2 /DNDEBUG" CACHE STRING "" FORCE)
    set(CMAKE_C_FLAGS_RELEASE "/O2 /DNDEBUG" CACHE STRING "" FORCE)
  endif()
endif()

if (CMAKE_C_COMPILER_ID MATCHES MSVC OR CMAKE_C_COMPILER_ID MATCHES Clang)
  if(CMAKE_BUILD_TYPE STREQUAL "Release")
    set(CMAKE_C_COMPILER_FLAGS "${CMAKE_C_COMPILER_ID} ${CMAKE_C_COMPILER} ${CMAKE_C_FLAGS_RELEASE} ${CMAKE_C_FLAGS}")
  else()
    set(CMAKE_C_COMPILER_FLAGS "${CMAKE_C_COMPILER_ID} ${CMAKE_C_COMPILER} ${CMAKE_C_FLAGS_DEBUG} ${CMAKE_C_FLAGS}")
  endif()
endif()

if(ENABLE_AMD_OPT)
  add_definitions(-DAMD_OPT_ALL)
endif()

if(ENABLE_AMD_MPIFFT)
  add_definitions(-DAMD_OPT_MPIFFT)
endif()

if(ENABLE_AMD_MPI_VADER_LIMIT)
  add_definitions(-DAMD_MPI_VADER_LIMIT_SET)
endif()

if(ENABLE_MPI)
  if(ENABLE_QUAD_PRECISION)
    message (FATAL_ERROR "quad precision is not supported in MPI")
  endif()
  find_package(MPI REQUIRED)
  set(HAVE_MPI TRUE)
endif()

if(ENABLE_AMD_TRANS)
    if(ENABLE_AMD_OPT)
  	add_definitions(-DAMD_OPT_TRANS)
    else()
	message(FATAL_ERROR "Main optimization switch ENABLE_AMD_OPT must be enabled to enable this option.")
    endif()
endif()

if (ENABLE_AMD_APP_OPT)
  if (NOT (ENABLE_MPI OR  ENABLE_QUAD_PRECISION OR ENABLE_LONG_DOUBLE))
    if(ENABLE_AMD_OPT)
      add_definitions(-DAMD_APP_OPT_LAYER)
    else ()
      message(FATAL_ERROR "Main optimization switch ENABLE_AMD_OPT must be enabled to enable this option.")
    endif ()
  else ()
    message(FATAL_ERROR "AMD application optimization layer is not supported for MPI execution and in Quad or Long double precisions.")
  endif ()
endif ()

if(ENABLE_AMD_DYNAMIC_DISPATCHER)
	message (FATAL_ERROR "AMD Dynamic Dispatcher feature is not available on windows.")
endif()

option (DISABLE_FORTRAN "Disable Fortran wrapper routines" OFF)

if (CMAKE_C_COMPILER_ID MATCHES MSVC OR CMAKE_C_COMPILER_ID MATCHES Clang)
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /MP")
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /W0 ")
endif ()

#print configurations
message("---cmake configurations---")
if (CMAKE_C_COMPILER_ID MATCHES Clang)
  message("AMD_ARCH selected:${AMD_ARCH}")
endif ()
message(CMAKE_C_COMPILER_ID : ${CMAKE_C_COMPILER_ID})
message(CMAKE_BUILD_TYPE : ${CMAKE_BUILD_TYPE})
message(BUILD_SHARED_LIBS : ${BUILD_SHARED_LIBS})
message(BUILD_TESTS : ${BUILD_TESTS})
message(ENABLE_OPENMP : ${ENABLE_OPENMP})
message(ENABLE_THREADS : ${ENABLE_THREADS})
message(WITH_COMBINED_THREADS : ${WITH_COMBINED_THREADS})
message(ENABLE_FLOAT : ${ENABLE_FLOAT})
message(ENABLE_LONG_DOUBLE : ${ENABLE_LONG_DOUBLE})
message(ENABLE_QUAD_PRECISION : ${ENABLE_QUAD_PRECISION})
message(ENABLE_SSE : ${ENABLE_SSE})
message(ENABLE_SSE2 : ${ENABLE_SSE2})
message(ENABLE_AVX : ${ENABLE_AVX})
message(ENABLE_AVX2 : ${ENABLE_AVX2})
message(ENABLE_AVX512 : ${ENABLE_AVX512})
message(ENABLE_AMD_OPT : ${ENABLE_AMD_OPT})
message(ENABLE_AMD_MPIFFT : ${ENABLE_AMD_MPIFFT})
message(ENABLE_AMD_MPI_VADER_LIMIT : ${ENABLE_AMD_MPI_VADER_LIMIT})
message(ENABLE_MPI : ${ENABLE_MPI})
message(ENABLE_AMD_TRANS : ${ENABLE_AMD_TRANS})
message(DISABLE_FORTRAN : ${DISABLE_FORTRAN})
message(CMAKE_CXX_FLAGS : ${CMAKE_CXX_FLAGS})
message(CMAKE_C_FLAGS : ${CMAKE_C_FLAGS})
message(CMAKE_C_FLAGS_RELEASE : ${CMAKE_C_FLAGS_RELEASE})
message(CMAKE_C_COMPILER_FLAGS : ${CMAKE_C_COMPILER_FLAGS})
message(ENABLE_AMD_FAST_PLANNER : ${ENABLE_AMD_FAST_PLANNER})
message(ENABLE_AMD_TOP_N_planner : ${ENABLE_AMD_TOP_N_planner})
message(ENABLE_AMD_APP_OPT : ${ENABLE_AMD_APP_OPT})

include(GNUInstallDirs)

include (CheckIncludeFile)
check_include_file (alloca.h         HAVE_ALLOCA_H)
check_include_file (altivec.h        HAVE_ALTIVEC_H)
check_include_file (c_asm.h          HAVE_C_ASM_H)
check_include_file (dlfcn.h          HAVE_DLFCN_H)
check_include_file (intrinsics.h     HAVE_INTRINSICS_H)
check_include_file (inttypes.h       HAVE_INTTYPES_H)
check_include_file (libintl.h        HAVE_LIBINTL_H)
check_include_file (limits.h         HAVE_LIMITS_H)
check_include_file (mach/mach_time.h HAVE_MACH_MACH_TIME_H)
check_include_file (malloc.h         HAVE_MALLOC_H)
check_include_file (memory.h         HAVE_MEMORY_H)
check_include_file (stddef.h         HAVE_STDDEF_H)
check_include_file (stdint.h         HAVE_STDINT_H)
check_include_file (stdlib.h         HAVE_STDLIB_H)
check_include_file (string.h         HAVE_STRING_H)
check_include_file (strings.h        HAVE_STRINGS_H)
check_include_file (sys/types.h      HAVE_SYS_TYPES_H)
check_include_file (sys/time.h       HAVE_SYS_TIME_H)
check_include_file (sys/stat.h       HAVE_SYS_STAT_H)
check_include_file (sys/sysctl.h     HAVE_SYS_SYSCTL_H)
check_include_file (time.h           HAVE_TIME_H)
check_include_file (uintptr.h        HAVE_UINTPTR_H)
check_include_file (unistd.h         HAVE_UNISTD_H)
if (HAVE_TIME_H AND HAVE_SYS_TIME_H)
  set (TIME_WITH_SYS_TIME TRUE)
endif ()

include (CheckPrototypeDefinition)
check_prototype_definition (drand48 "double drand48 (void)" "0" stdlib.h HAVE_DECL_DRAND48)
check_prototype_definition (srand48 "void srand48(long int seedval)" "0" stdlib.h HAVE_DECL_SRAND48)
check_prototype_definition (cosl "long double cosl( long double arg )" "0" math.h HAVE_DECL_COSL)
check_prototype_definition (sinl "long double sinl( long double arg )" "0" math.h HAVE_DECL_SINL)
check_prototype_definition (memalign "void *memalign(size_t alignment, size_t size)" "0" malloc.h HAVE_DECL_MEMALIGN)
check_prototype_definition (posix_memalign "int posix_memalign(void **memptr, size_t alignment, size_t size)" "0" stdlib.h HAVE_DECL_POSIX_MEMALIGN)

include (CheckSymbolExists)
check_symbol_exists (clock_gettime time.h HAVE_CLOCK_GETTIME)
check_symbol_exists (gettimeofday sys/time.h HAVE_GETTIMEOFDAY)
check_symbol_exists (getpagesize unistd.h HAVE_GETPAGESIZE)
check_symbol_exists (drand48 stdlib.h HAVE_DRAND48)
check_symbol_exists (srand48 stdlib.h HAVE_SRAND48)
check_symbol_exists (memalign malloc.h HAVE_MEMALIGN)
check_symbol_exists (posix_memalign stdlib.h HAVE_POSIX_MEMALIGN)
check_symbol_exists (mach_absolute_time mach/mach_time.h HAVE_MACH_ABSOLUTE_TIME)
check_symbol_exists (alloca alloca.h HAVE_ALLOCA)
if (NOT HAVE_ALLOCA)
  unset (HAVE_ALLOCA CACHE)
  check_symbol_exists (alloca malloc.h HAVE_ALLOCA)
endif ()
check_symbol_exists (isnan math.h HAVE_ISNAN)
check_symbol_exists (snprintf stdio.h HAVE_SNPRINTF)
check_symbol_exists (strchr string.h HAVE_STRCHR)
check_symbol_exists (sysctl unistd.h HAVE_SYSCTL)

if (UNIX)
  set (CMAKE_REQUIRED_LIBRARIES m)
endif ()
check_symbol_exists (cosl math.h HAVE_COSL)
check_symbol_exists (sinl math.h HAVE_SINL)

include (CheckTypeSize)
check_type_size ("float" SIZEOF_FLOAT)
check_type_size ("double" SIZEOF_DOUBLE)
check_type_size ("int" SIZEOF_INT)
check_type_size ("long" SIZEOF_LONG)
check_type_size ("long long" SIZEOF_LONG_LONG)
check_type_size ("unsigned int" SIZEOF_UNSIGNED_INT)
check_type_size ("unsigned long" SIZEOF_UNSIGNED_LONG)
check_type_size ("unsigned long long" SIZEOF_UNSIGNED_LONG_LONG)
check_type_size ("size_t" SIZEOF_SIZE_T)
check_type_size ("ptrdiff_t" SIZEOF_PTRDIFF_T)
math (EXPR SIZEOF_INT_BITS "8 * ${SIZEOF_INT}")
set (C_FFTW_R2R_KIND "C_INT${SIZEOF_INT_BITS}_T")

if (MSVC)
  add_definitions(-D_CRT_SECURE_NO_WARNINGS)
endif(MSVC)

string(TIMESTAMP TODAY "%Y%m%d")
add_compile_definitions(AOCL_FFTW_VERSION="AOCL-FFTW 4.2.1 Build ${TODAY}")

find_library (LIBM_LIBRARY NAMES m)
if (LIBM_LIBRARY)
  set (HAVE_LIBM TRUE)
endif ()


if (ENABLE_THREADS)
  find_package (Threads)
endif ()
if (Threads_FOUND)
  if(CMAKE_USE_PTHREADS_INIT)
    set (USING_POSIX_THREADS 1)
  endif ()
  set (HAVE_THREADS TRUE)
endif ()

if (WIN32)
    set(OpenMP_libomp_LIBRARY "C:/Program Files/LLVM/lib/libomp.lib" CACHE STRING "openmp library path")
endif (WIN32)

if (ENABLE_OPENMP)
  find_package (OpenMP)
endif ()
if (OPENMP_FOUND)
  set (HAVE_OPENMP TRUE)
endif ()

include (CheckCCompilerFlag)

if(CMAKE_C_COMPILER_ID MATCHES MSVC)
  if (ENABLE_AVX512)
    set (AVX512_FLAG /arch:AVX512)
    set(HAVE_AVX512 TRUE)
    message(STATUS "Enabling AVX512 instructions")
  endif()

  if (ENABLE_AVX2)
    set (AVX2_FLAG /arch:AVX2)
    set(HAVE_AVX2 TRUE)
    message(STATUS "Enabling AVX2 instructions")
    #not necessary to enable FMA in MSVC as its enabled
    #by default when AVX2 is enabled for x64
    # set (FMA_FLAG /arch:FMA)
    # add_definitions(/arch:FMA)
    # set(HAVE_FMA TRUE)
    # message(STATUS "Enabling FMA instructions")
  endif()

  if(ENABLE_AVX)
    set (AVX_FLAG /arch:AVX)
    set(HAVE_AVX TRUE)
    message(STATUS "Enabling AVX instructions")
  endif()

  if(ENABLE_SSE2)
    set (SSE2_FLAG /arch:SSE2)
    set(HAVE_SSE2 TRUE)
    message(STATUS "Enabling SSE2 instructions")
  endif()

  if(ENABLE_SSE)
    set (SSE_FLAG /arch:SSE)
    set(HAVE_SSE TRUE)
    message(STATUS "Enabling SSE instructions")
  endif()

  add_compile_definitions(${AVX_FLAG} ${AVX2_FLAG} ${SSE_FLAG} ${SSE2_FLAG})

elseif(CMAKE_C_COMPILER_ID MATCHES Clang)

  if (ENABLE_AVX512)
    set (AVX512_FLAG -mavx512f)
    set(HAVE_AVX512 TRUE)
    message(STATUS "Enabling AVX512 instructions")
  endif()

  if (ENABLE_AVX2)
    set (AVX2_FLAG -mavx2)
    set(HAVE_AVX2 TRUE)
    set (FMA_FLAG -mfma)
    set(HAVE_FMA TRUE)
    message(STATUS "Enabling AVX2 and FMA instructions")
  endif()

  if(ENABLE_AVX)
    set (AVX_FLAG -mavx)
    set(HAVE_AVX TRUE)
    message(STATUS "Enabling AVX instructions")
  endif()

  if(ENABLE_SSE2)
    set (SSE2_FLAG -msse2)
    set(HAVE_SSE2 TRUE)
    message(STATUS "Enabling SSE2 instructions")
  endif()

  if(ENABLE_SSE)
    set (SSE_FLAG -msse)
    set(HAVE_SSE TRUE)
    message(STATUS "Enabling SSE instructions")
  endif()

  add_compile_options(${AVX_FLAG} ${AVX2_FLAG} ${SSE_FLAG} ${SSE2_FLAG} ${FMA_FLAG})

else()

  if (ENABLE_SSE)
    foreach (FLAG "-msse" "/arch:SSE")
      unset (HAVE_SSE CACHE)
      unset (HAVE_SSE)
      check_c_compiler_flag (${FLAG} HAVE_SSE)
      if (HAVE_SSE)
        set (SSE_FLAG ${FLAG})
        break()
      endif ()
    endforeach ()
  endif ()

  if (ENABLE_SSE2)
    foreach (FLAG "-msse2" "/arch:SSE2")
      unset (HAVE_SSE2 CACHE)
      unset (HAVE_SSE2)
      check_c_compiler_flag (${FLAG} HAVE_SSE2)
      if (HAVE_SSE2)
        set (SSE2_FLAG ${FLAG})
        break()
      endif ()
    endforeach ()
  endif ()

  if (ENABLE_AVX)
    foreach (FLAG "-mavx" "/arch:AVX")
      unset (HAVE_AVX CACHE)
      unset (HAVE_AVX)
      check_c_compiler_flag (${FLAG} HAVE_AVX)
      if (HAVE_AVX)
        set (AVX_FLAG ${FLAG})
        break()
      endif ()
    endforeach ()
  endif ()

  if (ENABLE_AVX2)
    foreach (FLAG "-mavx2" "/arch:AVX2")
      unset (HAVE_AVX2 CACHE)
      unset (HAVE_AVX2)
      check_c_compiler_flag (${FLAG} HAVE_AVX2)
      if (HAVE_AVX2)
        set (AVX2_FLAG ${FLAG})
        break()
      endif ()
    endforeach ()
  endif ()

  # AVX2 codelets require FMA support as well
  if (ENABLE_AVX2)
    foreach (FLAG "-mfma" "/arch:FMA")
      unset (HAVE_FMA CACHE)
      unset (HAVE_FMA)
      check_c_compiler_flag (${FLAG} HAVE_FMA)
      if (HAVE_FMA)
        set (FMA_FLAG ${FLAG})
        break()
      endif ()
    endforeach ()
  endif ()

  if (ENABLE_AVX512)
    foreach (FLAG "-mavx512f" "/arch:AVX512")
      unset (HAVE_AVX512 CACHE)
      unset (HAVE_AVX512)
      check_c_compiler_flag (${FLAG} HAVE_AVX512)
      if (HAVE_AVX512)
        set (AVX512_FLAG ${FLAG})
        break()
      endif ()
    endforeach ()
  endif ()

endif ()

if (HAVE_SSE2 OR HAVE_AVX)
  set (HAVE_SIMD TRUE)
endif ()
file(GLOB           fftw_api_SOURCE                 api/*.c             api/*.h)
file(GLOB           fftw_dft_SOURCE                 dft/*.c             dft/*.h)
file(GLOB           fftw_dft_scalar_SOURCE          dft/scalar/*.c      dft/scalar/*.h)
file(GLOB           fftw_dft_scalar_codelets_SOURCE dft/scalar/codelets/*.c     dft/scalar/codelets/*.h)
file(GLOB           fftw_dft_simd_SOURCE            dft/simd/*.c        dft/simd/*.h)

file(GLOB           fftw_dft_simd_sse2_SOURCE       dft/simd/sse2/*.c   dft/simd/sse2/*.h)
file(GLOB           fftw_dft_simd_avx_SOURCE        dft/simd/avx/*.c    dft/simd/avx/*.h)
file(GLOB           fftw_dft_simd_avx2_SOURCE       dft/simd/avx2/*.c   dft/simd/avx2/*.h dft/simd/avx2-128/*.c   dft/simd/avx2-128/*.h)
file(GLOB           fftw_dft_simd_avx512_SOURCE     dft/simd/avx512/*.c dft/simd/avx512/*.h)
file(GLOB           fftw_kernel_SOURCE              kernel/*.c          kernel/*.h)
file(GLOB           fftw_rdft_SOURCE                rdft/*.c            rdft/*.h)
file(GLOB           fftw_rdft_scalar_SOURCE         rdft/scalar/*.c     rdft/scalar/*.h)

file(GLOB           fftw_rdft_scalar_r2cb_SOURCE    rdft/scalar/r2cb/*.c
                                                    rdft/scalar/r2cb/*.h)
file(GLOB           fftw_rdft_scalar_r2cf_SOURCE    rdft/scalar/r2cf/*.c
                                                    rdft/scalar/r2cf/*.h)
file(GLOB           fftw_rdft_scalar_r2r_SOURCE     rdft/scalar/r2r/*.c
                                                    rdft/scalar/r2r/*.h)

file(GLOB           fftw_rdft_simd_SOURCE           rdft/simd/*.c       rdft/simd/*.h)
file(GLOB           fftw_rdft_simd_sse2_SOURCE      rdft/simd/sse2/*.c  rdft/simd/sse2/*.h)
file(GLOB           fftw_rdft_simd_avx_SOURCE       rdft/simd/avx/*.c   rdft/simd/avx/*.h)
file(GLOB           fftw_rdft_simd_avx2_SOURCE      rdft/simd/avx2/*.c  rdft/simd/avx2/*.h rdft/simd/avx2-128/*.c  rdft/simd/avx2-128/*.h)
file(GLOB           fftw_rdft_simd_avx512_SOURCE    rdft/simd/avx512/*.c rdft/simd/avx512/*.h)

file(GLOB           fftw_reodft_SOURCE              reodft/*.c          reodft/*.h)
file(GLOB           fftw_simd_support_SOURCE        simd-support/*.c    simd-support/*.h)
file(GLOB           fftw_libbench2_SOURCE           libbench2/*.c       libbench2/*.h)

file(GLOB           fftw_mpi_SOURCE                 mpi/*.c             mpi/*.h)

list (REMOVE_ITEM   fftw_libbench2_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/libbench2/useropt.c)

list (REMOVE_ITEM   fftw_mpi_SOURCE ${CMAKE_SOURCE_DIR}/mpi/mpi-bench.c)
list (REMOVE_ITEM   fftw_mpi_SOURCE ${CMAKE_SOURCE_DIR}/mpi/testsched.c)

set(SOURCEFILES
    ${fftw_api_SOURCE}
    ${fftw_dft_SOURCE}
    ${fftw_dft_scalar_SOURCE}
    ${fftw_dft_scalar_codelets_SOURCE}
    ${fftw_dft_simd_SOURCE}
    ${fftw_kernel_SOURCE}
    ${fftw_rdft_SOURCE}
    ${fftw_rdft_scalar_SOURCE}

    ${fftw_rdft_scalar_r2cb_SOURCE}
    ${fftw_rdft_scalar_r2cf_SOURCE}
    ${fftw_rdft_scalar_r2r_SOURCE}

    ${fftw_rdft_simd_SOURCE}
    ${fftw_reodft_SOURCE}
    ${fftw_simd_support_SOURCE}
    ${fftw_threads_SOURCE}
)

set(fftw_par_SOURCE
    threads/api.c
    threads/conf.c
    threads/ct.c
    threads/dft-vrank-geq1.c
    threads/f77api.c
    threads/hc2hc.c
    threads/rdft-vrank-geq1.c
    threads/vrank-geq1-rdft2.c)

set (fftw_threads_SOURCE ${fftw_par_SOURCE} threads/threads.c)
set (fftw_omp_SOURCE ${fftw_par_SOURCE} threads/openmp.c)


include_directories (.)


if (WITH_COMBINED_THREADS)
  list (APPEND SOURCEFILES ${fftw_threads_SOURCE})
endif ()


if (HAVE_SSE2)
  list (APPEND SOURCEFILES ${fftw_dft_simd_sse2_SOURCE} ${fftw_rdft_simd_sse2_SOURCE})
endif ()

if (HAVE_AVX)
  list (APPEND SOURCEFILES ${fftw_dft_simd_avx_SOURCE} ${fftw_rdft_simd_avx_SOURCE})
endif ()

if (HAVE_AVX2)
  list (APPEND SOURCEFILES ${fftw_dft_simd_avx2_SOURCE} ${fftw_rdft_simd_avx2_SOURCE})
endif ()

if (HAVE_AVX512)
  set_source_files_properties(${fftw_dft_simd_avx512_SOURCE} ${fftw_rdft_simd_avx512_SOURCE} PROPERTIES COMPILE_OPTIONS "${AVX512_FLAG}")
  list (APPEND SOURCEFILES ${fftw_dft_simd_avx512_SOURCE} ${fftw_rdft_simd_avx512_SOURCE})
endif ()

set (FFTW_VERSION 3.3.10)

set (PREC_SUFFIX)
if (ENABLE_FLOAT)
  set (FFTW_SINGLE TRUE)
  set (BENCHFFT_SINGLE TRUE)
  set (PREC_SUFFIX f)
endif ()

if (ENABLE_LONG_DOUBLE)
  set (FFTW_LDOUBLE TRUE)
  set (BENCHFFT_LDOUBLE TRUE)
  set (PREC_SUFFIX l)
endif ()

if (ENABLE_QUAD_PRECISION)
  set (FFTW_QUAD TRUE)
  set (BENCHFFT_QUAD TRUE)
  set (PREC_SUFFIX q)
endif ()

if (ENABLE_AMD_FAST_PLANNER)
  if (NOT (ENABLE_QUAD_PRECISION OR ENABLE_LONG_DOUBLE))
    add_definitions(-DAMD_OPT_FAST_PLANNER)
  else ()
    message(FATAL_ERROR "AMD_FAST_PLANNER cannot be set for Quad and Long Double precision")
  endif ()
endif ()

if (ENABLE_AMD_TOP_N_planner)
  # Check if amd-top-n-planner is enabled with Quad and Long Double precision
  if (NOT (ENABLE_QUAD_PRECISION OR ENABLE_LONG_DOUBLE))
  # Check if both amd-fast-planner and amd-top-n-planner are enabled together
	if (NOT (ENABLE_AMD_FAST_PLANNER))
		# Check if amd-top-n-planner is enabled with mpi, openmp or threads
		if (NOT (ENABLE_MPI OR ENABLE_OPENMP OR ENABLE_THREADS))
			add_definitions(-DAMD_OPT_TOP_N_PLANNER)
		else ()
			message(FATAL_ERROR "AMD_TOP_N_PLANNER can not be enabled with mpi, openmp or threads as it is supported only for single threaded mode")
		endif ()
	else ()
    message(FATAL_ERROR "AMD_FAST_PLANNER and AMD_TOP_N_PLANNER can not be enabled together")
	endif ()
  else ()
    message(FATAL_ERROR "AMD_TOP_N_PLANNER cannot be set for Quad and Long Double precision")
  endif ()
  
  endif ()

set (fftw3_lib libfftw3${PREC_SUFFIX})

configure_file (cmake.config.h.in config.h @ONLY)
include_directories (${CMAKE_CURRENT_BINARY_DIR})

if (BUILD_SHARED_LIBS)
  add_definitions (-DFFTW_DLL)
  set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
endif ()

set(lib_name ${fftw3_lib})
add_library (${lib_name} ${SOURCEFILES})
target_include_directories (${lib_name} INTERFACE $<INSTALL_INTERFACE:include>)

if (Threads_FOUND)
  if (WITH_COMBINED_THREADS)
    set (lib_name ${fftw3_lib})
    add_library (${fftw3_lib} ${SOURCEFILES})
    target_include_directories (${lib_name} INTERFACE $<INSTALL_INTERFACE:include>)
    target_link_libraries (${lib_name} ${CMAKE_THREAD_LIBS_INIT})
  else ()
    set (lib_name ${fftw3_lib}_threads)
    add_library (${lib_name} ${fftw_threads_SOURCE} ${SOURCEFILES})
    target_include_directories (${lib_name} INTERFACE $<INSTALL_INTERFACE:include>)
    target_link_libraries (${lib_name} ${CMAKE_THREAD_LIBS_INIT})
  endif ()
endif ()
  
if (OPENMP_FOUND)
  set (lib_name ${fftw3_lib}_omp)
  add_library (${lib_name} ${fftw_omp_SOURCE} ${SOURCEFILES})
  target_include_directories (${lib_name} INTERFACE $<INSTALL_INTERFACE:include>)
  target_link_libraries (${lib_name} ${CMAKE_THREAD_LIBS_INIT})
  target_compile_options (${lib_name} PRIVATE ${OpenMP_C_FLAGS})  
  if (BUILD_SHARED_LIBS)
       target_link_libraries(${lib_name} PRIVATE OpenMP::OpenMP_C)
  endif(BUILD_SHARED_LIBS)
endif ()

if (HAVE_MPI)
  list (APPEND SOURCEFILES ${fftw_mpi_SOURCE})
  set (lib_name ${fftw3_lib}_mpi)
  add_library (${lib_name} ${fftw_omp_SOURCE} ${SOURCEFILES})
  target_include_directories (${lib_name} INTERFACE $<INSTALL_INTERFACE:include>)
  target_link_libraries (${lib_name} ${CMAKE_THREAD_LIBS_INIT})
  target_compile_options (${lib_name} PRIVATE ${OpenMP_C_FLAGS}) 
  if (BUILD_SHARED_LIBS AND OPENMP_FOUND)
       target_link_libraries(${lib_name} PRIVATE OpenMP::OpenMP_C)
  endif() 
endif ()


target_include_directories(${lib_name} PRIVATE ${CMAKE_SOURCE_DIR}/api)

if (CMAKE_C_COMPILER_ID MATCHES MSVC)
  target_compile_definitions (${lib_name} PRIVATE /bigobj)
endif ()

if (HAVE_SSE)
  target_compile_options (${lib_name} PRIVATE ${SSE_FLAG})
endif ()
if (HAVE_SSE2)
  target_compile_options (${lib_name} PRIVATE ${SSE2_FLAG})
endif ()
if (HAVE_AVX)
  target_compile_options (${lib_name} PRIVATE ${AVX_FLAG})
endif ()
if (HAVE_AVX2)
  target_compile_options (${lib_name} PRIVATE ${AVX2_FLAG})
endif ()
if (HAVE_FMA)
  target_compile_options (${lib_name} PRIVATE ${FMA_FLAG})
endif ()
#if (HAVE_AVX512)
#  target_compile_options (${lib_name} PRIVATE ${AVX512_FLAG})
#endif ()
if (HAVE_LIBM)
  target_link_libraries (${lib_name} m)
endif ()
if (HAVE_MPI)
  target_include_directories (${lib_name} PRIVATE ${MPI_INCLUDE_PATH})
  target_link_libraries (${lib_name} PRIVATE ${MPI_LIBRARIES})
endif ()

set_target_properties (${lib_name} PROPERTIES SOVERSION 3.5.7 VERSION 3)

install (TARGETS ${lib_name}
	  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
	  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
          ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

install(TARGETS ${lib_name}
          EXPORT FFTW3LibraryDepends
          RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
          LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
          ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

install (FILES api/fftw3.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

file(COPY api/fftw3.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/api/fftw3.f)
  install (FILES api/fftw3.f api/fftw3l.f03 api/fftw3q.f03 DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
  file(COPY api/fftw3.f api/fftw3l.f03 api/fftw3q.f03 DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
endif ()
if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/api/fftw3.f03.in)
  file (READ api/fftw3.f03.in FFTW3_F03_IN OFFSET 42)
  file (WRITE ${CMAKE_CURRENT_BINARY_DIR}/fftw3.f03 "! Generated automatically.  DO NOT EDIT!\n\n")
  file (APPEND ${CMAKE_CURRENT_BINARY_DIR}/fftw3.f03 "  integer, parameter :: C_FFTW_R2R_KIND = ${C_FFTW_R2R_KIND}\n\n")
  file (APPEND ${CMAKE_CURRENT_BINARY_DIR}/fftw3.f03 "${FFTW3_F03_IN}")
  install (FILES ${CMAKE_CURRENT_BINARY_DIR}/fftw3.f03 DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
  file(COPY ${CMAKE_CURRENT_BINARY_DIR}/fftw3.f03 DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
endif ()


if (HAVE_MPI)
	install (FILES mpi/fftw3-mpi.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
	if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/mpi/fftw3-mpi.f03.in)
	  file (READ mpi/fftw3-mpi.f03.in FFTW3-MPI_F03_IN OFFSET 42)
	  file (WRITE ${CMAKE_CURRENT_BINARY_DIR}/fftw3-mpi.f03 "! Generated automatically.  DO NOT EDIT!\n\n")
	  file (APPEND ${CMAKE_CURRENT_BINARY_DIR}/fftw3-mpi.f03 "${FFTW3-MPI_F03_IN}")
	  install (FILES ${CMAKE_CURRENT_BINARY_DIR}/fftw3-mpi.f03 DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
	endif ()
	if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/mpi/fftw3l-mpi.f03.in)
	  file (READ mpi/fftw3l-mpi.f03.in FFTW3L-MPI_F03_IN OFFSET 42)
	  file (WRITE ${CMAKE_CURRENT_BINARY_DIR}/fftw3l-mpi.f03 "! Generated automatically.  DO NOT EDIT!\n\n")
	  file (APPEND ${CMAKE_CURRENT_BINARY_DIR}/fftw3l-mpi.f03 "${FFTW3L-MPI_F03_IN}")
	  install (FILES ${CMAKE_CURRENT_BINARY_DIR}/fftw3l-mpi.f03 DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
	endif ()
	file(COPY mpi/fftw3-mpi.h ${CMAKE_CURRENT_BINARY_DIR}/fftw3-mpi.f03 ${CMAKE_CURRENT_BINARY_DIR}/fftw3l-mpi.f03 DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
endif ()

add_library (libbench2 STATIC ${fftw_libbench2_SOURCE})
add_library (libtestbench STATIC tests/bench.c tests/hook.c tests/fftw-bench.c tests/fftw-bench.h)

if(Threads_FOUND)
        set(CMAKE_C_COMPILER_CFLAGS ${PTHREAD_CFLAGS})
        if(NOT COMBINED_THREADS)
                set(LIBFFTWTHREADS libfftw3${PREC_SUFFIX}_threads)
        endif()
elseif(OPENMP_FOUND)
        set(CMAKE_C_COMPILER_CFLAGS ${OpenMP_C_FLAGS})
        set(LIBFFTWTHREADS libfftw3${PREC_SUFFIX}_omp)
endif()

set(fftw${PREC_SUFFIX}-wisdom_src tools/fftw-wisdom.c)

add_executable(fftw${PREC_SUFFIX}-wisdom ${fftw${PREC_SUFFIX}-wisdom_src})
target_link_libraries(fftw${PREC_SUFFIX}-wisdom libbench2  ${LIBFFTWTHREADS} ${lib_name} libtestbench)
if (OPENMP_FOUND)
    target_link_libraries(fftw${PREC_SUFFIX}-wisdom OpenMP::OpenMP_C)
endif ()
install(TARGETS fftw${PREC_SUFFIX}-wisdom ${INSTALL_TARGETS_DEFAULT_ARGS} DESTINATION ${CMAKE_INSTALL_BINDIR})


if (BUILD_TESTS)
  add_executable (bench tests/bench.c tests/hook.c tests/fftw-bench.c)
  target_link_libraries (bench libbench2)
  target_link_libraries (bench ${lib_name})
  if (OPENMP_FOUND)
    target_link_libraries(bench OpenMP::OpenMP_C)
  endif ()

  if (HAVE_MPI)
    add_executable (mpi-bench  mpi/mpi-bench.c tests/hook.c tests/fftw-bench.c)
    target_include_directories (mpi-bench PRIVATE ${CMAKE_SOURCE_DIR}/api)
    target_include_directories (mpi-bench PRIVATE ${MPI_INCLUDE_PATH})
    target_link_libraries (mpi-bench libbench2 ${MPI_LIBRARIES} ${lib_name})
    if (OPENMP_FOUND)
      target_link_libraries(mpi-bench OpenMP::OpenMP_C)
    endif ()
  endif ()

  enable_testing ()

  if (Threads_FOUND)

    macro (fftw_add_test problem)
      add_test (NAME ${problem} COMMAND bench -s ${problem})
    endmacro ()

    fftw_add_test (32x64)
    fftw_add_test (ib256)

  endif ()
endif ()

# pkgconfig file
set (prefix ${CMAKE_INSTALL_PREFIX})
set (exec_prefix ${CMAKE_INSTALL_PREFIX})
set (libdir ${CMAKE_INSTALL_FULL_LIBDIR})
set (includedir ${CMAKE_INSTALL_FULL_INCLUDEDIR})
set (VERSION ${FFTW_VERSION})
configure_file (fftw.pc.in fftw3${PREC_SUFFIX}.pc @ONLY)
install (FILES
          ${CMAKE_CURRENT_BINARY_DIR}/fftw3${PREC_SUFFIX}.pc
         DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig
         COMPONENT Development)

# cmake file
set (FFTW3_LIBRARIES "FFTW3::${lib_name}")
configure_file (FFTW3Config.cmake.in FFTW3${PREC_SUFFIX}Config.cmake @ONLY)
configure_file (FFTW3ConfigVersion.cmake.in FFTW3${PREC_SUFFIX}ConfigVersion.cmake @ONLY)
install (FILES
          ${CMAKE_CURRENT_BINARY_DIR}/FFTW3${PREC_SUFFIX}Config.cmake
          ${CMAKE_CURRENT_BINARY_DIR}/FFTW3${PREC_SUFFIX}ConfigVersion.cmake
	  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/fftw3${PREC_SUFFIX}
         COMPONENT Development)

export (TARGETS ${lib_name} NAMESPACE FFTW3:: FILE ${PROJECT_BINARY_DIR}/FFTW3LibraryDepends.cmake)
install(EXPORT FFTW3LibraryDepends
        NAMESPACE FFTW3::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/fftw3${PREC_SUFFIX}
        COMPONENT Development)
