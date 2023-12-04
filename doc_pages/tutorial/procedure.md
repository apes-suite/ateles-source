title: General Procedure for Running the Testcases

For running one of the testcases:

1.  Go to the corresponding testcase folder:
    -   `ateles/tutorial/testcase/euler` for the Euler testcase, or
    -   `ateles/tutorial/testcase/maxwell` for the Maxwell testcase.
1.  In this folder there is already a lua input file `ateles.lua`
    in which the configuration for the simulation is defined.
1.  Create a folder called `restart`
    (the input file contains a section for configuring the restart options
    including the loaction for storing the restart files
    which is set to `write = './restart/'` by default).
    If you do not create said folder, the simulation will fail
    and give the following error:
    `File Open MPI WARNING: MPI_ERR_NO_SUCH_FILE: no such file or directory`.
1.  Run the simulation with `../../../build/ateles ateles.lua`,
    where `../../../build/ateles` is the path to the executable
    and `ateles.lua` is the name of the input file.

Now, you can see the simluation output in the terminal
and finally the data files will be stored in the restart folder.
