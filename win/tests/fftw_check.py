import os
import sys
import subprocess
import yaml


class FftwCheck:

    @staticmethod
    def check_execution():
        """
        :Method Name: check_execution
        :Description: reads the command from input file and executes if the respective executables are present
        :parameter  : None
        :return     : None
        """
        try:
            with open(r'commands.yaml') as file:
                input_file = yaml.safe_load(file) \

                try:
                    if (sys.argv[1] == '') or (sys.argv[1] == "--h") or (sys.argv[1] == "--help"):
                        print("Below options are available \n\nUse '*_mt' options for Multithreaded"
                              " build validation\nUse '*_mpi'"
                              "options for MPI build validation\n")
                        print("usage: python fftw_check.py ",
                              end='[check | smallcheck | bigcheck | ')
                        for var in input_file.keys():
                            print(var, end=' | ')
                        print(' --h | --help]')
                        sys.exit()
                except IndexError:
                    print("Below options are available \n\n"
                          "Use '*_mt' options for Multithreaded build validation\nUse "
                          "'*_mpi' options for MPI build validation\n")
                    print("usage: python fftw_check.py ", end='[check | smallcheck | bigcheck | ')
                    for var in input_file.keys():
                        print(var, end=' | ')
                    print(' --h | --help]')
                    sys.exit()

                try:
                    mpi_run = False
                    if os.path.exists("bench.exe"):
                        if os.path.exists("mpi-bench.exe"):
                            mpi_run = True
                            numcheck = 10 if sys.argv[1] == "check" \
                                else 100 if sys.argv[1] == "bigcheck" else 2
                            check_list = ['_st', '_mt', '_mpi'] \
                                if mpi_run else ['_st', '_mt']
                        else:
                            check_list = ['_st', '_mt'] \
                                if mpi_run else ['_st', '_mt']
                    elif not 'mpi' in sys.argv[1]:
                        print("bench.exe is not present ")
                        sys.exit()
                    elif os.path.exists("mpi-bench.exe"):
                        mpi_run = True
                        numcheck = 10 if sys.argv[1] == "check" \
                            else 100 if sys.argv[1] == "bigcheck" else 2
                    else:
                        print("bench.exe is not present ")
                        sys.exit()

                    temp = 'basic' if sys.argv[1] == "check" \
                        else 'big' if sys.argv[1] == "bigcheck" else 'a few'
                    if (sys.argv[1] == "check") or (sys.argv[1] == "smallcheck") \
                            or (sys.argv[1] == "bigcheck"):

                        for command in range(len(check_list)):
                            print("=" * 50, sys.argv[1] + check_list[command],
                                  "execution", "=" * 50)
                            for i in range(len(input_file[sys.argv[1] + check_list[command]])):
                                process = subprocess.Popen(input_file[sys.argv[1]
                                                                      + check_list[command]][i],
                                                           bufsize=1,
                                                           universal_newlines=True,
                                                           stdout=subprocess.PIPE,
                                                           stderr=subprocess.STDOUT)
                                for line in iter(process.stdout.readline, ''):
                                    print(line[:-1])

                                    sys.stdout.flush()
                                process.wait()
                                # errcode = process.returncode

                                if check_list[command] == '_mpi':
                                    if '--nthreads' in input_file[sys.argv[1]
                                                                  + check_list[command]][i]:
                                        print("-" * 80, "\n", "\t" * 2,
                                              "MPI FFTW threaded transforms passed "
                                              "{} tests!".format(numcheck),
                                              "\n", "-" * 80)
                                    else:
                                        print("-" * 80, "\n", "\t" * 2,
                                              "MPI FFTW transforms passed "
                                              "{} tests, {} CPU".format(numcheck, i + 1),
                                              "\n", "-" * 80)

                            if check_list[command] == '_st':
                                print("*" * 80, "\n", "\t" * 2, "FFTW "
                                                                "transforms passed %s "
                                                                "tests!\n" % temp, "*" * 80)
                            elif check_list[command] == '_mt':
                                print("*" * 80, "\n", "\t" * 2, "FFTW "
                                                                "threaded transforms passed %s "
                                                                "tests!\n" % temp,
                                      "*" * 80)
                    else:
                        if 'mpi' in sys.argv[1]:
                            if mpi_run:
                                pass
                            else:
                                print('mpi-bench.exe is not present '
                                      ', please select other valid arguments')
                                sys.exit()

                        print("=" * 50, sys.argv[1], "execution", "=" * 50)
                        for i in range(len(input_file[sys.argv[1]])):
                            process = subprocess.Popen(input_file[sys.argv[1]][i], bufsize=1,
                                                       universal_newlines=True,
                                                       stdout=subprocess.PIPE,
                                                       stderr=subprocess.STDOUT)

                            for line in iter(process.stdout.readline, ''):
                                print(line[:-1])
                                sys.stdout.flush()
                            process.wait()
                            # errcode = process.returncode

                            if '_mpi' in sys.argv[1]:
                                if "--nthreads" in input_file[sys.argv[1]][i]:
                                    print("-" * 80, "\n", "\t" * 2,
                                          "MPI FFTW threaded transforms "
                                          "passed {} tests!".format(numcheck), "\n",
                                          "-" * 80)
                                else:
                                    print("-" * 80, "\n", "\t" * 2,
                                          "MPI FFTW transforms passed "
                                          "{} tests, {} CPU".format(numcheck, i + 1), "\n",
                                          "-" * 80)

                        temp = 'a few' if ("smallcheck" in sys.argv[1]) else 'big' if (
                                "bigcheck" in sys.argv[1]) else 'basic'
                        if '_st' in sys.argv[1]:
                            print("*" * 80, "\n", "\t" * 2, "FFTW "
                                                            "transforms passed %s "
                                                            "tests!\n" % temp, "*" * 80)
                        elif '_mt' in sys.argv[1]:
                            print("*" * 80, "\n", "\t" * 2, "FFTW "
                                                            "threaded transforms passed %s "
                                                            "tests!\n" % temp,
                                  "*" * 80)
                except UnboundLocalError as error:
                    print(error)

        except KeyError:
            print("\nPlease enter a valid argument")
            print("usage: python fftw_check.py ", end='[check | smallcheck | bigcheck |')
            for var in input_file.keys():
                print(var, end=' | ')
            print(' --h | --help]')


if __name__ == "__main__":
    FftwCheck.check_execution()
