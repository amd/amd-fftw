AMD OPTIMIZED FFTW
------------------

AMD Optimized FFTW is the optimized FFTW implementation targeted for 
AMD EPYC CPUs. It is developed on top of FFTW (version fftw-3.3.8).
All known features and functionalities of FFTW are retained and supported
as it is with this AMD optimized FFTW library.

AMD Optimized FFTW achieves higher performance than the FFTW 3.3.8 due to
various optimizations involving improved SIMD Kernel functions, improved copy
functions (cpy2d and cpy2d_pair used in rank-0 transform and buffering plan),
improved 256-bit kernels selection by Planner and an optional in-place 
transpose for large problem sizes. AMD Optimized FFTW improves the performance
of in-place MPI FFT over FFTW 3.3.8 by employing a faster in-place MPI
transpose function. AMD Optimized FFTW provides a new fast planner as an
extension to the original planner that improves planning time of various
planning modes in general and PATIENT mode in particular. As of AMD FFTW 3.0.1,
a new feature called Top N planner is introduced that minimizes single-threaded
run-to-run variations.

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

INSTALLATION FROM AMD Optimized FFTW GIT REPOSITORY:

After downloading the latest stable release from the git repository,
https://github.com/amd/amd-fftw, follow the below steps to configure and
build it for AMD EPYC processor based on Naples, Rome and future 
generation architectures.

     ./configure --enable-sse2 --enable-avx --enable-avx2 
                 --enable-mpi --enable-openmp --enable-shared 
                 --enable-amd-opt --enable-amd-mpifft 
                 --prefix=<your-install-dir>
     make
     make install

The configure option "--enable-amd-opt" enables all the improvements and 
optimizations targeted for AMD EPYC CPUs.
When enabling configure option "--enable-amd-opt", do not use the 
configure option "--enable-generic-simd128" or "--enable-generic-simd256".

The optional configure option "--enable-amd-mpifft" enables the MPI FFT
related optimizations.

An optional configure option "--enable-amd-mpi-vader-limit" is supported that 
controls enabling of AMD's new MPI transpose algorithms. When using this 
configure option, the user needs to set --mca btl_vader_eager_limit
appropriately (current preference is 65536) in the MPIRUN command.

The new fast planner can be enabled using optional configure option 
"--enable-amd-fast-planner". It is supported for single and double precisions.

An optional configure option "AMD_ARCH" is supported that can be set to CPU 
architecture values like "auto" or "znver1" or "znver2" for AMD EPYC processors.

An optional configure option "--enable-amd-trans" is provided that may benefit
the performance of transpose operations in case of very large FFT problem sizes.
This is by default not enabled and provided as an experimental optional switch. 

By default, configure script enables double-precision mode. User should pass
appropriate configure options to enable the single-precision or quad-precision
or long-double mode.

CONTACTS
--------

AMD Optimized FFTW is developed and maintained by AMD.
You can contact us on the email-id aoclsupport@amd.com.
You can also raise any issue/suggestion on the git-hub repository at
https://github.com/amd/amd-fftw/issues

ACKNOWLEDGEMENTS
----------------

FFTW was developed by Matteo Frigo and Steven G. Johnson. We thank Matteo Frigo
for his support provided to us.
