#FFTW check execution script

check execution script covers check, smallcheck, bigcheck, paranoid-check and exhaustive-check

##Requirements
* Install latest version of python from python.org(preferably python 3.5 or greater)
* Add python path and scripts path to the environment variable path
* Install python PyYAML module using the following command
  pip install PyYAML
*Install the latest version of perl from the below giver URL
  http://strawberryperl.com/
*Add perl path to the environment variable path

#Copy all the files present in win/tests directory to the directory where .exe is present
#Open the command prompt and execute the python script and provide an argument(check_st or check_mt or smallcheck_mpi etc)
#To know about the valid arguments , use help option(python fftw_check.py --help)
  For example:
    python fftw_check.py smallcheck
    python fftw_check.py check
    python fftw_check.py bigcheck_mt
    python fftw_check.py smallcheck_st

#Output can be seen on the command prompt


#FFTW benchmark execution script
For help on executing .bat file run below command from command prompt
fftw_st_benchmark.bat --help
fftw_mt_benchmark.bat --help
fftw_mpi_benchmark.bat --help

Sample commands:
fftw_st_benchmark.bat --verify bench.exe
fftw_st_benchmark.bat --verify-only bench.exe
fftw_mt_benchmark.bat --verify bench.exe
fftw_mt_benchmark.bat --verify-only bench.exe
fftw_mpi_benchmark.bat --verify-only mpi-bench.exe
