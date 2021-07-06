@echo OFF
setlocal EnableDelayedExpansion
set NL=^


REM two empty line required
REM FFTW-MPI Lib Batch Test
REM This shell script benchmarks a program on a series of test problems
REM
REM #set -ix
set me=%0


if "%1" == "" (
   echo Try '!me! --help' for more information.
   exit /b )

set usage_options=--h --help --verify --verify-only --verify-tolerance --accuracy --accuracy-rounds --k --keep-going --time-min --max --maxnd -o --user-option

set usage_txt=   Usage: !me! [OPTION] program !NL! ^
  -h, --help         print this help, then exit !NL! ^
  --verify           verify each transform before timing !NL! ^
  --verify-only      verify each transform but do not time !NL! ^
  --verify-tolerance set error tolerance for --verify !NL! ^
  --accuracy         run accuracy test !NL! ^
  --accuracy-rounds  set number of rounds for --accuracy !NL! ^
  -k, --keep-going   continue after verification error !NL! ^
  --time-min:X       set minimum measurement time !NL! ^
  --maxn:N           set maximum allowed problem size !NL! ^
  --maxnd:N          set maximum allowed multi-dimensional problem size !NL! ^
  -o:X, --user-option:X    undocumented !NL!

set verify=no
set accuracy=no
set useropt=""
set speed=yes
set maxnd=1073741824 #1 billion size
set maxn=67777216
set time_min=""
set keep_going=no
set tolerance=""
set rounds=""
set arounds=""

:GetLastArg
set /a argCount+=1
set program=%~1
shift
if not "%~1"=="" goto GetLastArg

set delim_char=":"
set arg_cnt=0

for %%i in (%*) do (
   set /a arg_cnt = !arg_cnt!+1
   set optargset=yes
   set arg_val=%%~i

   if !arg_val! == !program! goto NEXT
   echo %%~i|find %delim_char% >nul
   if errorlevel 1 set optargset=no
   for /F "tokens=1,2 delims=: " %%a in ("%%i") do (
      set vararg=%%~a
      set optarg=%%~b
      )
   REM echo vararg=!vararg!
   if !optargset! == yes (
       echo !usage_options!|find "!vararg!" >nul
       if errorlevel 1 (
          echo Invalid option !vararg!
          exit /b )

       if !vararg! == --verify-rounds (
       set rounds=!optarg!
       )
       if !vararg! == --verify-tolerance (
       set tolerance=!optarg!
       )
       if !vararg! == --keep-going (
       set keep_going=!optarg!
       )
       if !vararg! == -k (
       set keep_going=!optarg!
       )
       if !vararg! == --accuracy-rounds (
       set arounds=!optarg!
       )
       if !vararg! == --time-min (
       set time_min=!optarg!
       )
       if !vararg! == --maxn (
       set maxn=!optarg!
       )
       if !vararg! == --maxnd (
       set maxnd=!optarg!
       )
       if !vararg! == --user-option (
       set useropt="%useropt% --user-option=!optarg!"
       )
    ) else (
       if !vararg! == --verify-rounds (
       echo error: missing argument to !vararg!
       exit /b
       )
       if !vararg! == --verify-tolerance (
       echo error: missing argument to !vararg!
       exit /b
       )
       if !vararg! == --accuracy-rounds (
       echo error: missing argument to !vararg!
       exit /b
       )
       if !vararg! == --time-min (
       echo error: missing argument to !vararg!
       exit /b
       )
       if !vararg! == --maxn (
       echo error: missing argument to !vararg!
       exit /b
       )
       if !vararg! == --maxnd (
       echo error: missing argument to !vararg!
       exit /b
       )
       if !vararg! == --keep-going (
       echo error: missing argument to !vararg!
       exit /b
       )
       if !vararg! == -k (
       echo error: missing argument to !vararg!
       exit /b
       )
       if !vararg! == --user-option (
       echo error: missing argument to !vararg!
       exit /b
       )
       if !arg_val! == --help (
       echo !usage_txt!
       exit /b
       )
       if !arg_val! == --h (
       echo !usage_txt!
       exit /b
       )
       if !arg_val! == --keep-going (
       set keep-going=yes
       )
       if !arg_val! == --k (
       set keep-going=yes
       )
       if !arg_val! == --accuracy (
       set accuracy=yes
       set speed=no
       )
       if !arg_val! == --verify (
       set verify=yes
       )
       if !arg_val! == --verify-only (
       set verify=yes
       set speed=no
       )
       echo !usage_options!|find "!arg_val!" >nul
       if errorlevel 1 (
          echo Invalid option !vararg!
          exit /b )
   )
)

:NEXT
if !arg_cnt! == 1 (
   if !arg_val! == --help echo !usage_txt!
   if !arg_val! == -h echo !usage_txt!
   if !arg_val! == !program! goto :NEXT_1
   echo Try '!me! --help' for more information.
   exit /b )

:NEXT_1
REM MPIFFT 1D Sizes
set SIZES_1D=524288 1048576 2097152 4194304 8388608 16777216 33554432 67108864 390625 823543 4782969 7962624 1000000 8000000 16000000 32000000 64000000

REM MPIFFT 2D Sizes
set SIZES_2D=2048x2048 2048x4096 4096x4096 8192x8192 8192x16384 729x2187 625x3125 6561x15625 27000x27000 1000x10000 10000x10000

REM MPIFFT 3D Sizes
set SIZES_3D=256x256x256 512x512x512 512x512x1024 1024x1024x512 1024x1024x1024 243x243x243 343x343x343 243x243x729 625x625x625 200x180x216 300x1000x100 1000x1000x1000

set DIRECTIONS=f
set ALL_SIZES=!SIZES_1D! !SIZES_2D! !SIZES_3D!
REM set ALL_SIZES=2 4x4 1960x1960 2048x2048 4x4x4 16x1024x64

set PLACE=i
set REALITY=c
set CORES=16 32 64 128

if not exist "!program!" (
   echo !program! file does not exist
   exit /b
)

REM test "$speed" = "no" || test -n "$time_min" || time_min=`$program --print-time-min`
if !speed! == no (
for /f %%i in ('%program% --print-time-min') do set time_min=%%~i )

if !time_min! == "" (
   for /f %%i in ('!program! --print-time-min') do set time_min=%%~i
)

REM precision=`$program --print-precision`
for /f %%i in ('!program! --print-precision') do set precision=%%~i

REM shorten the name
if !precision! == single set precision=s
if !precision! == double set precision=d

REM name=`$program --info name`
for /f %%i in ('!program! --info name') do set name=%%~i

if not %tolerance%=="" (
set tolerance=--verify-tolerance %tolerance%
set vflags=!vflags! !tolerance!
)

if not !rounds! == "" (
set rounds=--verify-rounds !rounds!
set vflags=!vflags! !rounds!
)

if not !arounds! == "" (
   set arounds=--accuracy-rounds !arounds!
   set vflags=!vflags! !arounds!
   set aflags=!arounds!
)

for %%r in (%CORES%) do (
set rank_var=%%~r
for %%p in (%PLACE%) do (
   set plc_var=%%~p
   for %%r in (%REALITY%) do (
      set r_val=%%~r
      for %%s in (%ALL_SIZES%) do (
         REM :NEXT_LOOP
         set size_var=%%~s
         set size_val=%%~s
         set oned=no
         REM #for core in $CORES; do
         set curmaxn=%maxnd%
         echo !size_var!|find "x">nul
         if errorlevel 1 set curmaxn=%maxn%

         echo !size_var!|find "x">nul
         if errorlevel 1 set oned=yes

         set size_val=!size_val:x=*!
         set /a size_val=!size_val!

         if !size_val! LSS !curmaxn! (
             for %%d in (%DIRECTIONS%) do (
                set d_val=%%~d
                set problem=!plc_var!!r_val!!d_val!!size_var!
                REM #doable=`$program $useropt --can-do $problem`
                REM set doable="#t"
                set print=yes
                set /A ITERS=3
                for /l %%v in (1,1,!ITERS!) do (
                   set acc=""
                   if !verify! == yes (
                      for /f "tokens=*" %%i in ('!program! !vflags! --verbose --verify !problem!') do set acc=%%i
                      REM echo !acc!
                      if !acc! == 'FAILED FAILED FAILED' if not !keep_going! == yes exit /b
                    ) else ( set acc="" )
                   if !speed! == yes (
                      if !useropt! == "" set useropt=
                      if !rank_var! == 16 (
                         for /f "tokens=*" %%j in ('mpiexec --map-by L3cache --bind-to core -np !rank_var! !program! !useropt! --report-benchmark --time-min !time_min! -opatient -r500 --speed !problem!') do set time_val=%%j
                         if !time_val! == "" set time_val="FAILED FAILED"
                      )
                      if !rank_var! == 32 (
                         for /f "tokens=*" %%j in ('mpiexec --map-by L3cache --bind-to core -np !rank_var!  !program! !useropt! --report-benchmark --time-min !time_min! -opatient -r500 --speed !problem!') do set time_val=%%j
                         if !time_val! == "" set time_val="FAILED FAILED"
                      ) else (
                      if !rank_var! == 64 (
                         for /f "tokens=*" %%j in ('mpiexec --map-by core --bind-to core -np !rank_var! !program! !useropt! --report-benchmark --time-min !time_min! -opatient -r500 --speed !problem!') do set time_val=%%j
                         if !time_val! == "" set time_val="FAILED FAILED"
                      )
                      )
                    ) else ( set time_val="" )
                   if !accuracy! == yes (
                      if !oned! == yes (
                        for /f "tokens=*" %%k in ('!program! !useropt! !aflags! --accuracy !problem!') do set acc=%%k
                    ) else ( set print="no" )
                   )
                   if !print! == yes (
                     echo !rank_var! !name! !precision!!r_val!!plc_var!!d_val! !size_var! !time_val! !acc! )
                )
                )
           )
        )
    )
)
)
endlocal