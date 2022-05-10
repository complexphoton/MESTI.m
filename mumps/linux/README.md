## Prerequisites for MUMPS installation on a Linux

We need the following tools before compiling MUMPS and its MATLAB interface on the Linux.

### MATLAB 

The MATLAB should be installed before compiling MUMPS-MATLAB interface. The MATLAB compiler `mex` is required when compiling the MUMPS-MATLAB interface. 

If MATLAB is already installed on the Linux cluster with the Lmod module system, you can type  
```shell
module load matlab
```

If the Linux cluster does not use Lmod module system, you may add the path of your MATLAB to `PATH` by
```shell
export PATH=".../matlab/bin/:$PATH"
```
where `...` is the path to MATLAB on the cluster.

### GNU compilers collection

The compilation of MUMPS requires both C and Fortran compilers. Both C and Fortran compliers are included in GNU Compiler Collection (GCC), which should be installed on a cluster. On a local machine, you will need to install them by yourself.

On a cluster, it is likely that GCC has already been installed. You can type
```shell
gcc -v
```
to check if GCC is available.

If GCC is not loaded or a different version is preferred, you can load GCC on a Lmod module system by typing
```shell
module load gcc
```
In non-Lmod module system, you may type,
```shell
.../gcc/setup.sh
```
where `...` is the path to GCC on the cluster.

On a local machine with Debian-based distributions, you can run
```shell
sudo apt install build-essential
```
to install the C compiler and some other compiling tools, such as <code>make</code>. You will also need to run
```shell
sudo apt install gfortran
```
to install <code>gfortran</code>.

### BLAS and LAPACK

MUMPS requires both BLAS and LAPACK libraries, which are standard libraries on Linux cluster. These libraries are also included in many implementations, such as MKL and OpenBLAS. 

In the example below, we use MKL on the cluster and OpenBLAS on the local machine. MKL is faster than OpenBLAS for running MESTI on an Intel CPU, though it requires more storage space. If you would like to use MKL on a local machine, please download the installer [here](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.html). Make sure to use custom installation and choose MKL only. The default installation takes 7.5 GB while MKL takes 1.5 GB (for version 2022.2.0).

Different clusters install MKL in different paths. You will need the correct path for the linker to find corresponding BLAS and LAPACK libraries. In the provided `Makefile.inc`, we assume the MKL path has been exported to an environment variable called `MKLROOT`. Here we show how to export the correct `MKLROOT`. The BLAS and LAPACK libraries can be found under `$(MKLROOT)/lib/intel64`. 

If you are using the Lmod module system and MKL is installed, you can use 

```shell
module load intel
module load intel-mkl
echo $MKLROOT
```

`MKLROOT`  is assigned and should be printed out, if MKL is available. The `MKLROOT` should be under the path `.../mkl`, where `...` is the Intel directory path.

If Lmod module system is not used in your cluster, then you can also try to find the path like `.../mkl`, where `...` is the Intel directory path. Then you can assign the `MKLROOT` by 

```shell
source .../mkl/bin/mklvars.sh intel64
echo $MKLROOT
```

`MKLROOT`  is assigned and should be printed out.

In some cases, the cluster cannot find the BLAS and LAPACK libraries by itself when you run MATLAB interface for MUMPS. You will also need to override MATLAB's own MKL. To solve those issues, you can append those library paths to `LD_PRELOAD` and `LD_LIBRARY_PATH` before running MATLAB. For example, if BLAS and LAPACK libraries are used through MKL, you can type,

```shell
export LD_PRELOAD=$LD_PRELOAD:$MKLROOT/lib/intel64/libmkl_intel_lp64.so:$MKLROOT/lib/intel64/libmkl_sequential.so:$MKLROOT/lib/intel64/libmkl_intel_thread.so:$MKLROOT/lib/intel64/libmkl_core.so
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MKLROOT/lib/intel64
```

If you still cannot figure out the path of the libraries, you may need to contact the maintenance team of your cluster.

For OpenBLAS, you can first install it by running

```shell
sudo apt install libopenblas-dev
```

in the terminal and then use the corresponding `LIBBLAS` and `LAPACK` in the `Makefile.inc` and `make.inc` we provided.