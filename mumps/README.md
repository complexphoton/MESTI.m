# MUMPS Installation

The SCSA method in MESTI requires the MUMPS package. Here are steps to install it.

## Download MUMPS
Please go to the [MUMPS website](http://mumps.enseeiht.fr/index.php?page=dwnld) and fill out the form to request the MUMPS source code.

## Compile MUMPS
After obtaining the MUMPS source code, read the <code>INSTALL</code> file to choose and modify the <code>Makefile.inc</code> to fit your environment and machine. You'll need to modify the environmental variables in <code>Makefile.inc</code> such as the paths of the BLAS library (required), LAPACK library (required), and METIS ordering library (optional): <code>LIBBLAS</code>, <code>LAPACK</code>, and <code>LMETIS</code>.

Two examples of <code>Makefile.inc</code> we use ourselves are provided below, one for installation on a Linux cluster, one for installation on MacOS. Note you still need to modify the environmental variables based on your machine.

1. [Linux cluster](./linux/Makefile.inc)
2. [MacOS](./macOS/Makefile.inc)

In you are on a cluster, you may need to load the modules for the libraries used.

Then, type

<code>make d z</code>

which will compile sequential MUMPS with double precision for real and complex variables (*i.e.*, <code>dmumps</code> and <code>zmumps</code>). 

If there is no error, go to the <code>lib</code> folder to check if the following files have been generated: <code>libdmumps.a</code> (or <code>libdmumps.so</code>) and <code>libzmumps.a</code> (or <code>libzmumps.so</code>).

## Compile MATLAB interface for MUMPS

Go to the <code>MATLAB</code> folder inside the MUMPS folder and modify <code>make.inc</code>.  Similarly, you'll need to modify the environmental variables such as the paths of the MUMPS folder (required), BLAS library (required), and the METIS ordering library (optional): <code>MUMPS_DIR</code>, <code>LIBBLAS</code>, and <code>LMETIS</code>.

Two examples of <code>make.inc</code> are provided below:
1. [Linux cluster](./linux/make.inc)
2. [MacOS](./macOS/make.inc)

When simulating large systems, there can be segmentation fault due to bugs in the MATLAB interface. We suggest replacing the original <code>mumpsmex.c</code> with this modified one, [mumpsmex.c](mumpsmex.c), which modifies four lines such that it disables reading the scaling array from MATLAB, and it no longer outputs the scaling array to MATLAB; these lines are where the segmentation fault happens. Since MESTI does not use user-specified scaling arrays, these modifications do not affect functionality.

In you are on a cluster, you may need to load the software module for MATLAB.

Then, type

<code>make</code>

which will compile the MATLAB interface for MUMPS. Check that files <code>dmumpsmex.mexa64</code> and <code>zmumpsmex.mexa64</code> have been generated in the <code>MATLAB</code> folder.

You can now add this folder to the search path using the <code>addpath</code> command in MATLAB. 


## Test MUMPS
You can run the examples in the MATLAB interface folder of MUMPS: <code>simple_example.m</code>, <code>multiplerhs_example.m</code>, <code>sparserhs_example.m</code>, <code>schur_example.m</code>, and <code>zsimple_example.m</code>. If any of them does not run successfully, check the compilation of MUMPS or the MATLAB interface.
