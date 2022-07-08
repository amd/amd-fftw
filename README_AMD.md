AOCL-FFTW
---------

AOCL-FFTW is AMD optimized version of FFTW implementation targeted for 
AMD EPYC CPUs. It is developed on top of FFTW (version fftw-3.3.10).
All known features and functionalities of FFTW are retained and supported
as it is with this AMD optimized FFTW library.

AOCL-FFTW achieves higher performance than the original FFTW 3.3.10 due to its
various optimizations involving improved SIMD Kernel functions, improved copy
functions (cpy2d and cpy2d_pair used in rank-0 transform and buffering plan),
improved 256-bit kernels selection by Planner and an optional in-place 
transpose for large problem sizes. AOCL-FFTW improves the performance
of in-place MPI FFTs over FFTW 3.3.10 by employing a faster in-place MPI
transpose function. AOCL-FFTW provides a new fast planner mode as an
extension to the original planner that improves planning time of various
planning modes in general and PATIENT mode in particular. Another new planning
mode called Top N planner is also available that minimizes single-threaded
run-to-run variations. AOCL-FFTW has a feature called AMD's application
optimization layer that speeds up HPC and scientific applications. As of 
AOCL-FFTW-3.2, a new feature called dynamic dispatcher is implemented that 
can build a single portable optimized library for execution on a wide range of
x86 CPU architectures.

FFTW is a free collection of fast C routines for computing the
Discrete Fourier Transform and various special cases thereof in one or more
dimensions. It includes complex, real, symmetric, and parallel transforms, 
and can handle arbitrary array sizes efficiently.

The doc/ directory contains the manual in texinfo, PDF, info, and HTML
formats.  Frequently asked questions and answers can be found in the
doc/FAQ/ directory in ASCII and HTML.

For a quick introduction to calling FFTW, see the "Tutorial" section
of the manual.

INSTALLATION
------------

INSTALLATION FROM AOCL-FFTW GIT REPOSITORY:

After downloading the latest stable release from the git repository,
https://github.com/amd/amd-fftw, follow the below steps to configure and
build it for AMD EPYC processor based on Naples, Rome, Milan and future 
generation architectures.

     ./configure --enable-sse2 --enable-avx --enable-avx2 
                 --enable-mpi --enable-openmp --enable-shared 
                 --enable-amd-opt --enable-amd-mpifft 
                 --prefix=<your-install-dir>
     make
     make install

The configure option "--enable-amd-opt" enables all the improvements and 
optimizations targeted for AMD EPYC CPUs. For enabling various optional
configure options provided for AMD EPYC CPUs, the master optimization switch
"--enable-amd-opt" must be kept enabled.

When enabling configure option "--enable-amd-opt", do not use the 
configure option "--enable-generic-simd128" or "--enable-generic-simd256".

The optional configure option "--enable-amd-mpifft" enables the MPI FFT
related optimizations.

An optional configure option "--enable-amd-mpi-vader-limit" is supported that 
controls enabling of AMD's new MPI transpose algorithms. When using this 
configure option, the user needs to set --mca btl_vader_eager_limit
appropriately (current preference is 65536) in the MPIRUN command.

The new fast planner can be enabled using optional configure option 
"--enable-amd-fast-planner". It is supported in single and double precisions.

Top N planner mode can be enabled using optional configure option
"--enable-amd-top-n-planner" to minimize run-to-run variations in performance.
It is supported in single-threaded execution in single and double precisions.

An optional configure option "AMD_ARCH" is supported that can be set to CPU 
architecture values like "auto" or "znver1" or "znver2" or "znver3" for AMD 
EPYC processors.

The optional configure option "--enable-amd-app-opt" turns on AMD's application
optimization layer to benefit performance of HPC and scientific applications.
Currently it is developed for complex and real (r2c and c2r) DFT problem types
in double and single precisions. It is not supported for MPI FFTs, r2r real DFT
problem types, Quad or Long double precisions, and split array format.

Dynamic dispatcher achieves Function Multi-versioning by using compiler's
attributes. Use "--enable-dynamic-dispatcher" configure option to enable this
feature. It is supported for GCC compiler and Linux based systems for now.
The set of x86 CPUs on which the single portable library can work depends upon
the highest level of CPU SIMD instruction set with which it is configured.

An optional configure option "--enable-amd-trans" is provided that may benefit
the performance of transpose operations in case of very large FFT problem sizes.
This is by default not enabled and provided as an experimental optional switch. 

By default, configure script enables double-precision mode. User should pass
appropriate configure options to enable the single-precision or quad-precision
or long-double mode.

CONTACTS
--------

AOCL-FFTW is developed and maintained by AMD.
You can contact us on the email-id aoclsupport@amd.com.
You can also raise any issue/suggestion on the git-hub repository at
https://github.com/amd/amd-fftw/issues

ACKNOWLEDGEMENTS
----------------

FFTW was developed by Matteo Frigo and Steven G. Johnson. We thank Matteo Frigo
for his support provided to us.
