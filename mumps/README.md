# MUMPS

Some of our methods can use the MUMPS solver. Here we provide steps to download and install MUMPS.

## Download MUMPS
Please go to the [MUMPS website](http://mumps-solver.org/) and fill out the form to get the MUMPS source code. After then, you should receive the MUMPS source code.

## Compile MUMPS
After getting the MUMPS source code, please go through the <code>INSTALL</code> file to choose and modify the <code>Makefile.inc</code> to fit your environment and machine.

The rule of thumb is that the users should modify the environmental variables in <code>Makefile.inc</code> to compile MUMPS successfully. For example, the path of BLAS library (required), LAPACK library (required), and METIS library (optional): <code>LIBBLAS</code>, <code>LAPACK</code>, and <code>LMETIS</code>.

Two of the examples of <code>Makefile.inc</code> is provided below for Linux cluster and MacOS systems, but users may still need to modify the environmental variables according to the local machine.

1. [Linux cluster](./linux/Makefile.inc)
2. [MacOS](./macOS/Makefile.inc)

After properly modifying the <code>Makefile.inc</code>, and then type

<code>make d z</code>

The process is to compile sequential MUMPS with double precision and double complex versions (i.e. dmumps and zmumps). 

If there is no error, then go to <code>lib</code> directory to check if there are  <code>libdmumps.a</code> (or <code>libdmumps.so</code>) and <code>libzmumps.a</code> (or <code>libzmumps.so</code>). Then it is fine to go to next step.

## Compile MATLAB interface to MUMPS

Go to <code>MATLAB</code> directory inside the <code>MUMPS_5.4.1</code> directory and modify <code>make.inc</code>. 

Similarly, the rule of thumb is that the users should modify the environmental variables in <code>make.inc</code> to compile MATLAB interface to MUMPS successfully. For example, the path of MUMPS directory (required), BLAS library (required), and METIS library (optional): <code>MUMPS_DIR</code>, <code>LIBBLAS</code>, and <code>LMETIS</code>.


Two of the examples of <code>make.inc</code> is provided below for Linux cluster and MacOS systems, but users may still need to modify the environmental variables according to the local machine.

1. [Linux cluster](./linux/make.inc)
2. [MacOS](./macOS/make.inc)

For large systems the MATLAB interface could cause segmentation fault. We suggest replacing the original <code>mumpsmex.c</code> with the modified file here. The new <code>mumpsmex.c</code> modifies four lines such that it disables reading scaling array from the user and it no longer outputs the scaling array to MATLAB. This appears to be the lines where the segmentation fault happens. MESTI do not use user-specified scaling arrays thus those modifications will not affect its functionality. Please refer to MUMPS user guide if you plan to use MUMPS for other packages. The modified <code>mumpsmex.c</code> can be found in the following link. 

[mumpsmex.c](mumpsmex.c)

After modifying make.inc and replacing <code>mumpsmex.c</code>, preload MATLAB, and then type

<code>make</code>

The process is to compile the MATLAB interface to MUMPS. Check whether there are <code>dmumpsmex.mexa64</code>, <code>zmumpsmex.mexa64</code> in <code>MATLAB</code> directory.


## Test MATLAB interface to MUMPS
Preload all the libraries needed in MUMPS. Users can run the five examples in <code>MUMPS_5.4.1/MATLAB</code> directory: <code>simple_example.m</code>, <code>multiplerhs_example.m</code>, <code>sparserhs_example.m</code>, <code>schur_example.m</code>, and <code>zsimple_example.m</code> to check that everything runs smoothly. If one of them cannot run successfully, please check the compilation for MUMPS or MATLAB interface.


## Preload MATLAB interface to MUMPS
If everything works well, then it is ready to call MUMPS solver in mesti(). The users just need to <code>addpath YOUR_DIRECTORY/MUMPS_5.4.1/MATLAB</code> on the head of the MATLAB script. Note that <code>YOUR_DIRECTORY/MUMPS_5.4.1/MATLAB</code> is the path for MATLAB interface to MUMPS.
