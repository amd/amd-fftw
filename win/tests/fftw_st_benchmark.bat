@echo OFF
setlocal EnableDelayedExpansion
set NL=^


REM two empty line required
REM FFTW-ST Lib Batch Test
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
set maxn=16777216
set maxnd=1048576
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
set SIZES_1D=2 4 6 8 9 12 15 16 18 24 32 36 64 80 108 128 210 256 504 512 1000 1024 1960 2048 4096 4725 8192 10368 16384 27000 32768 65536 75600 131072 165375 262144 362880 524288 1048576 1594323 1953125 2097152 4194304 4782969 5764801 8388608

set SIZES_2D=4x4 5x5 6x6 7x7 8x8 4x8 8x4 9x9 10x10 11x11 12x12 13x13 14x14 15x15 16x16 25x24 32x32 48x48 49x49 60x60 72x56 64x64 75x75 80x80 84x84 128x64 16x512 96x96 100x100 105x105 112x112 120x120 128x128 144x144 180x180 512x64 256x128 240x240 256x256 64x1024 360x360 512x512 1000x1000 1024x1024 1960x1960 2048x2048 3360x3360 4096x4096 4725x4725 8192x8192 10368x10368 16384x16384 27000x27000 32768x32768

set SIZES_3D=4x4x4 5x5x5 6x6x6 7x7x7 8x8x8 9x9x9 10x10x10 11x11x11 12x12x12 13x13x13 14x14x14 15x15x15 16x16x16 4x8x16 24x25x28 32x32x32 48x48x48 49x49x49 60x60x60 72x60x56 64x64x64 75x75x75 80x80x80 256x64x32 84x84x84 96x96x96 100x100x100 16x1024x64 105x105x105 112x112x112 120x120x120 128x128x128 144x144x144 512x128x64 180x180x180 256x128x256 240x240x240 256x256x256 512x64x1024 360x360x360 512x512x512

set DIRECTIONS=f
set ALL_SIZES=!SIZES_1D! !SIZES_2D! !SIZES_3D!

set PLACE=i
set REALITY=c
set CORES=2 4 8 16 32 64 128

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
                      for /f "tokens=*" %%f in ('!program! !useropt! --report-benchmark --time-min !time_min! -opatient --speed !problem!') do set time_val=%%f
                      if !time_val! == 'FAILED FAILED' if not !keep_going! == yes exit /b
                    ) else ( set time_val="" )
                   if !accuracy! == yes (
                      if !oned! == yes (
                        for /f "tokens=*" %%k in ('!program! !useropt! !aflags! --accuracy !problem!') do set acc=%%k
                    ) else ( set print="no" )
                   )
                   if !print! == yes echo !name! !precision!!r_val!!plc_var!!d_val! !size_var! !time_val! !acc!
                )
                )
           )
        )
    )
)
endlocal