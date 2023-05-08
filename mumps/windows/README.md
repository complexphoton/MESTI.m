## Prerequisites for MUMPS installation on Windows

We need the following tools before compiling MUMPS and its MATLAB interface on Windows. 

### Windows Subsystem for Linux (WSL)

Directly building MUMPS on Windows can be very challenging. To overcome this issue, it is highly recommended to use Windows Subsystem for Linux (WSL) as a viable solution. 

Please install WSL following the instruction [here](https://learn.microsoft.com/en-us/windows/wsl/install). 

### MATLAB on WSL

The MATLAB should be installed before compiling MUMPS-MATLAB interface. The MATLAB compiler `mex` is required when compiling the MUMPS-MATLAB interface. 

Please follow the standard process to download and install MATLAB on WSL. Alternatively, MATLAB Package Manager ([mpm](https://github.com/mathworks-ref-arch/matlab-dockerfile/blob/main/MPM.md)) can be used to install MATLAB.

After MATLAB is installed on the WSL, you can add the path of your MATLAB to `PATH` by
```shell
export PATH=".../matlab/bin/:$PATH"
```
where `...` is the path to MATLAB on the WSL. If you use mpm to install MATLAB, the default path is  `/usr/share/matlab/bin/`.

### GNU compilers collection on WSL

The compilation of MUMPS requires both C and Fortran compilers. Both C and Fortran compilers are included in GNU Compiler Collection (GCC). On WSL, you will need to install them by yourself.

On WSL, you can run
```shell
sudo apt install build-essential
```
to install the C compiler and some other compiling tools, such as <code>make</code>. You will also need to run
```shell
sudo apt install gfortran
```
to install <code>gfortran</code>.

### BLAS and LAPACK on WSL

MUMPS requires both BLAS and LAPACK libraries. These libraries are also included in many implementations, such as MKL and OpenBLAS. 

In the example below, we use MKL on WSL. Please download the installer of MKL [here](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.html). Make sure to use custom installation and choose MKL only.

You will need the correct path for the linker to find corresponding BLAS and LAPACK libraries. In the provided `Makefile.inc`, we assume the MKL path has been exported to an environment variable called `MKLROOT`. Here we show how to export the correct `MKLROOT`. The BLAS and LAPACK libraries can be found under `$(MKLROOT)/lib/intel64`. 

You can assign the `MKLROOT` by (for version 2023.1.0)

```shell
source /opt/intel/oneapi/mkl/2023.1.0/env/vars.sh
echo $MKLROOT
```

`MKLROOT`  is assigned and should be printed out.

In some cases, the WSL cannot find the BLAS and LAPACK libraries by itself when you run MATLAB interface for MUMPS. You will also need to override MATLAB's own MKL. To solve those issues, you can append those library paths to `LD_PRELOAD` and `LD_LIBRARY_PATH` before running MATLAB. For example, if BLAS and LAPACK libraries are used through MKL, you can type,

```shell
export LD_PRELOAD=$LD_PRELOAD:$MKLROOT/lib/intel64/libmkl_intel_lp64.so:$MKLROOT/lib/intel64/libmkl_gnu_thread.so:$MKLROOT/lib/intel64/libmkl_core.so
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MKLROOT/lib/intel64
```

Note that MUMPS supports shared memory, multithreaded parallelism through the use of multithreaded
BLAS libraries.  We provide `Makefile.inc` and `make.inc` on windows and they can activate the OpenMP feature. You can use the environment variable `OMP_NUM_THREADS` to set the number of threads.

### Note on schur_example.m on WSL

When running the test script `schur_example.m`, MATLAB may abort. This error seems to MATLAB, not from MUMPS. This can be ignored. As long as we can run the other test scripts (`simple_example.m`, `multiplerhs_example.m`, `sparserhs_example.m`, and `zsimple_example.m`) smoothly, we should be able to run MESTI.m and the [examples](../../examples).