# MUMPS Installation
MESTI.m uses the sequential version of MUMPS for the augmented partial factorization (APF) method, and optionally for the factorize-and-solve method. Here are steps to install MUMPS.

## Download MUMPS
Go to the [MUMPS website](https://mumps-solver.org/index.php?page=dwnld) and fill out the download request form. The MUMPS maintainers will email you the download link.

## Prerequisites
To compile the sequential version of MUMPS and its MATLAB interface, you need compilation tools like <code>make</code>, <code>ar</code>, and <code>ranlib</code>, C and Fortran compilers, BLAS library, LAPACK library, and a compatible <code>mex</code> compiler. Instructions specific to the operating system are provided below:
1. [Linux](./linux)
2. [macOS](./macOS)
3. [Windows](./windows)

If memory usage is important for you, you can optionally install the [METIS](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview) (version 5.1.0) program for graph partitioning (not to be confused with MESTI). Download it (you can just download archive/metis-5.1.0.tar.gz), decompress it, and set it to double precision (use <code>#define REALTYPEWIDTH 64</code> in <code>include/metis.h</code>). Compile it with <code>make config; make</code>. If you have write access to <code>/usr/local</code>, you can install METIS to there with <code>sudo make install</code>; otherwise you can move the METIS folder to where you want and specify its path when compiling MUMPS. Later, if you set <code>opts.use_METIS = true</code> in <code>mesti()</code> or <code>mesti2s()</code>, MUMPS will use METIS instead of the default AMD method for matrix ordering. From our experience, in 2D, AMD is usually faster when using the APF method (because its analysis stage is faster), but METIS can sometimes reduce memory usage. Which one is better depends on the problem. In 3D, the factorization time generally dominates over the analysis time, and the METIS ordering is much better (both in memmory usage and in factorization time) than AMD.

## Compile MUMPS
Suppose you downloaded the 5.6.0 version of MUMPS to your ~/Downloads/ folder. Then, go to the folder where you want to compile MUMPS, and enter
```
tar zxvf ~/Downloads/MUMPS_5.6.0.tar.gz
cd MUMPS_5.6.0
```
in terminal.

Read the file <code>INSTALL</code>, copy the closest <code>Makefile.inc</code> file in the <code>Make.inc</code> folder, and modify it to fit your environment and machine. Most importantly, in <code>Makefile.inc</code> you need to specify:
 - <code>CC</code>: the C compiler to use
 - <code>FC</code>: the Fortran compiler to use
 - <code>FL</code>: the Fortran linker to use
 - <code>LAPACK</code>: how the Fortran compiler can link to the LAPACK library
 - <code>LIBBLAS</code>: how the Fortran compiler can link to the BLAS library

If you installed METIS, you also need to specify
 - <code>LMETISDIR</code>: path to the folder where the METIS library is
 - <code>IMETIS</code>: path to the folder where the METIS headers are

and add <code>-Dmetis</code> to <code>ORDERINGSF</code>.

Examples of <code>Makefile.inc</code> are provided below:
1. [Linux](./linux/Makefile.inc)
2. [macOS](./macOS/Makefile.inc)
3. [Windows](./windows/Makefile.inc)

To download, click the link above, and click on the "Download raw file" button to the right of the "Raw" button.

These <code>Makefile.inc</code> examples for Linux and Windows enable multithreading by linking to a multithreaded BLAS library and compiling MUMPS with OpenMP activated.

The <code>Makefile.inc</code> example for macOS does not use multithreading, because there we link to Apple's [vecLib](https://developer.apple.com/documentation/accelerate/veclib) within its Accelerate framework for the BLAS library. On Macs with an Apple Silicon processor, the Accelerate framework utilizes the Apple Matrix coprocessor (AMX), so it is not multithreaded but achieves speed comparable to or faster than multithreaded OpenBLAS (*e.g.*, see [this benchmark](https://github.com/OpenMathLib/OpenBLAS/issues/2814#issuecomment-771505972)). Multithreading is, therefore, not necessary.

If you have a Mac with an Intel processor, OpenBLAS would be faster than Accelerate. In that case, you can comment out the line <code>LIBBLAS = -framework Accelerate</code> in <code>Makefile.inc</code> by adding a <code>#</code> in the beginning of the line, and uncomment the  <code>#LIBBLAS = -L/usr/local/opt/openblas/lib</code> line by deleting its <code>#</code>. In this case, multithreading can provide further speed-up, but we did not work out a multithreaded <code>Makefile.inc</code> example since Apple has moved away from Intel processors.

After done with <code>Makefile.inc</code>, enter
```
make d z
```
in terminal, which will compile the sequential version of MUMPS with double precision for real and complex variables (*i.e.*, <code>dmumps</code> and <code>zmumps</code>).

If there is no error, check if the following files have been generated in the <code>lib</code> folder: <code>libdmumps.a</code> (or <code>libdmumps.so</code>) and <code>libzmumps.a</code> (or <code>libzmumps.so</code>).

Warning messages from the Fortran compiler are normal and can be ignored.

If there is error, read the message and try to figure out where it comes from and/or look it up and address it. Before recompiling with <code>make d z</code>, be sure to type <code>make clean</code> first to remove files generate by the previous compile attempt.

## Compile the MATLAB interface for MUMPS

<code>cd</code> into the <code>MATLAB</code> folder inside the MUMPS folder and modify <code>make.inc</code>.  Most importantly, you need to specify
 - <code>MEX</code>: the <code>mex</code> compiler to use
 - <code>MUMPS_DIR</code>: path to the folder where MUMPS is compiled
 - <code>LIBFORT</code>: how <code>mex</code> can link to the Fortran library
 - <code>LIBBLAS</code>: how <code>mex</code> can link to the BLAS and LAPACK libraries

If you installed METIS, you also need to specify
 - <code>LMETISDIR</code>: path to the folder where the METIS library is

and uncomment the <code>#LMETIS     = -L$(LMETISDIR) -lmetis</code> line.

Examples of <code>make.inc</code> are provided below:
1. [Linux](./linux/make.inc)
2. [macOS](./macOS/make.inc)
3. [Windows](./windows/make.inc)

**IMPORTANT**: When simulating large systems (such as the [metalens example](../examples/2d_metalens)), one may encounter segmentation fault due to a bug in the MATLAB interface that comes with MUMPS. Please replace the original <code>mumpsmex.c</code> in this <code>MATLAB</code> folder with a modified one here, [mumpsmex.c](mumpsmex.c), which modifies four lines to disable reading the scaling array from MATLAB and returning the scaling array to MATLAB; these lines are where the segmentation fault happens. Since MESTI does not use user-specified scaling arrays, these modifications do not affect functionality.

Then, enter
```
make
```
in terminal, which will compile the MATLAB interface for MUMPS. Check that files <code>dmumpsmex.mexa64</code> and <code>zmumpsmex.mexa64</code> (or <code>\*.mexmaca64</code> on a Mac) have been generated in the <code>MATLAB</code> folder.


## Test MUMPS

Now, open MATLAB, <code>cd</code> to the MATLAB interface folder of MUMPS, and run the following test scripts:
- <code>simple_example.m</code>
- <code>multiplerhs_example.m</code>
- <code>sparserhs_example.m</code>
- <code>schur_example.m</code>
- <code>zsimple_example.m</code>.

If any of them does not run successfully, look back at the compilation of MUMPS or the MATLAB interface to see if there were serious warning messages.

If they all pass, congratulations! You are done.

Add this folder of MATLAB interface for MUMPS to the search path of MATLAB using the <code>addpath</code> command in MATLAB, so <code>mesti()</code> can find the function file <code>zmumps.m</code> there. You can also have this path to be added every time MATLAB starts by editing <code>startup.m</code> with the <code>edit(fullfile(userpath,'startup.m'))</code> command.

If you want to use METIS but your cluster cannot find the METIS library when you run the MATLAB interface for MUMPS, you can append the METIS library to your <code>LD_LIBRARYP_PATH</code> by

```shell
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$LMETISDIR
```

`LMETISDIR` is the path to the folder where the METIS library is.
