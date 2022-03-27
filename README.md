# MESTI.m

## Overview
**MESTI.m**  is a finite-difference frequency-domain (FDFD) package implementing **Schur complement scattering analysis (SCSA)** algorithm, which is described in [arXiv2203.xxxxx](https://arxiv.org/abs/2203.xxxxx). MESTI.m (an acronym for Maxwell's Equations Solver with Thousands of Inputs) enables users to solve the Maxwell's equations with multiple inputs by orders-of-magnitude speedup and reduced memory usage.

MESTI.m is written in MATLAB and specifically solves the scalar wave equation in 2D transverse-magnetic wave. 3D vectorial MESTI version is under development in Julia language and will be available soon.  


MESTI.m includes two main interfaces for users to use the package: <code>mesti()</code> and <code>mesti2s()</code>.

- <code>mesti()</code> computes field profile or scattering matrix with multiple inputs and deals with the outgoing boundary conditions through perfectly matched layer (PML).

- <code>mesti2s()</code> computes field profile or the scattering matrix with multiple inputs in two-sided geometry and can deal with the outgoing boundary conditions through perfectly matched layer (PML) or self-energy.


## Installation

MESTI.m is written in MATLAB and is not required specifically installation as long as users have installed MATLAB 2019b or later in their machine. 

However, to fully utilize SCSA algorithm, users need to install MUMPS package. Currently, MESTI.m implement SCSA algorithm through calling MUMP as solver. Without MUMPS, MESTI.m would not really solve problems through SCSA.

For the better performance and fully utilize the feature of MESTI.m, install MUMPS is required and highly recommended. Please look at the following page for the steps to install MUMPS: 

[MUMPS installation](./mumps/README.md)


## Available Functions 
The functions in MESTI.m can be separated into two parts: external functions and internal functions. Users are expected to use the external functions as interfaces to run MESTI.m. On the other hand, users may not need to use the internal functions directly. 

### External functions
-   <code>mesti()</code> solves frequency-domain scattering problem in general 2D geometry.
-   <code>mesti2s()</code> solves frequency-domain scattering problem specifically in a two-sided geometry.
-   <code>mesti_build_channels()</code> sets up properties of channels in the homogeneous space.

### Internal functions
-   <code>mesti_build_fdfd_matrix.m</code> builds up the finite-difference frequency-domain operator in 2D.
-   <code>mesti_matrix_solver.m</code> computes matrices CA<sup>-1</sup>B or A<sup>-1</sup>B.

For more details about the functions, please see the documentations in the codes. Or users can see more details about in each function by typing help and the function name in MATLAB, for example,
<code>help mesti2s</code>. Users can also go through some examples to get familiar with the external functions.


## Examples

There are several examples available <code>examples/*</code>. They are listed in the following:

1D
1. [Fabry-Perot etalon](./examples/1d_fabry_perot/fabry_perot.md)
2. [Distributed Bragg reflector](./examples/1d_distributed_bragg_reflector/distributed_bragg_reflector.md)

2D 
1. [Open Channel Through Disorder](./examples/2d_open_channel_through_disorder/open_channel_through_disorder.md)
2. [Meta-atom in metalens](./examples/2d_metalens_meta_atom/metalens_meta_atom.md)

Users can go through the examples to get familiar with the package.


## Galleries
1. Wave propagation in Fabry-Perot etalon
<img src="./examples/1d_fabry_perot/fabry_perot_field_profile.gif" width="400" height="300"> 

2. Open channel through disorder
<img src="./examples/2d_open_channel_through_disorder/disorder_open_channel.gif" width="400" height="300"> 

3. Intensity profile for a mm-wide hyperbolic metalens 
<img src="./examples/2d_metalens_full/metalens_intensity_profile_0_degree.png" width="400" height="300"> 

4. Transmission efficiency and Strehl ratio scanning all angles for a mm-wide hyperbolic metalens 
<img src="./examples/2d_metalens_full/metalens_Strehl_ratio_and_transmission_eff.png" width="400" height="300"> 
