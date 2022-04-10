## Prerequisite 

You will need to install the following packages first for compiling MUMPS on macOS.

### Xcode command-line tools

The command-line tools have several utilities for compilation, such as <code>make</code> and <code>clang</code>. To install command-line tools type

```
xcode-select --install
```
in the terminal and follow the dialogs that open.

### Homebrew

Homebrew will be helpful for the installation of <code>gfortran</code> and <code>openblas</code> (optional) later. Homebrew can be installed by pasting the following line in the terminal.

```
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```

### GFortran

The compilation of MUMPS requires both C and Fortran compilers. The C compiler has been installed by the command-line tools. The Fortran compiler can be install by <code>brew</code>. After installing Homebrew, simply type 

```
brew install gfortran
```

in the terminal to install <code>gfortran</code>.

## Optional packages

The following packages are not required if you use the <code>Makefile.inc</code> and <code>make.inc</code> we provide. You can install them for comparing the performance of MUMPS with different BLAS implementations.

### OpenBLAS

We use the BLAS and LAPACK implementation from Apple's Accelerate framework by default, which does not require external installation. However, on an intel-based system it has been tested that OpenBLAS has better performance. Thus we also provide the option of linking with OpenBLAS in those make files. To install OpenBLAS, type

```
brew install openblas
```
in the terminal. And it will install OpenBLAS under the directory <code>/usr/local/opt/openblas</code> by default.