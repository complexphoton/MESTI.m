# MESTI

**MESTI** (Maxwell's Equations Solver with Thousands of Inputs) is a finite-difference frequency-domain (FDFD) solver that implements the **Schur complement scattering analysis (SCSA)** method described in [arXiv:2203.xxxxx](https://arxiv.org/abs/2203.xxxxx). It allows users to solve Maxwell's equations with thousands of inputs while using computing time and memory that is even less than what a typical direct solver uses to solve with just one input.

MESTI.m is written in MATLAB and considers transverse-magnetic waves in 2D. A 3D vectorial version of MESTI written in Julia is under development and will be released later.  

The user can specify arbitrary permittivity profiles of the system, arbitrary lists of input sources (user-specified or automatically built), and arbitrary lists of output projections (or no projection, in which case the complete field profiles are returned). MESTI implements all of the common boundary conditions, perfectly matched layer (PML) with both imaginary and real coordinate stretching, as well as exact outgoing boundaries in two-sided or one-sided geometries. In addition to SCSA, MESTI also implements conventional direct methods.

## Installation

No installation is required for MESTI itself; just download it and add the <code>MESTI.m/src</code> folder to the search path using the <code>addpath</code> command in MATLAB. The MATLAB version should be R2019b or later. (Using an earlier version is possible but requires minor edits.)

However, to use the SCSA method (which is the key distinguishing feature of MESTI), the user needs to install the MUMPS package and its MATLAB interface. Without MUMPS, MESTI will still run but will only use slower conventional solvers that are not based on the Schur complement. So, MUMPS installation is strongly recommended.  See this [MUMPS installation](./mumps) page for steps to install MUMPS.

## Functions 

The function <code>mesti(syst, B, C, D)</code> provides the most flexibility. The user can specify arbitrary permittivity profiles, arbitrary boundary conditions, PML on any or all sides, with any lists of sources and any lists of output projections (or no output projection) placed at any location. But the user needs to provide the source and/or projection profiles.

The function <code>mesti2s(syst, in, out)</code> deals specifically with two-sided or one-sided geometries where the boundary condition in y is closed (*e.g.*, periodic or PEC)  and the boundary condition in x is outgoing. The user only needs to specify which channels and/or wavefronts are of interest; <code>mesti2s()</code> builds the lists of input sources and output projection profiles, and then calls <code>mesti()</code> for the computation. <code>mesti2s()</code> also offers the additional features of (1) exact outgoing boundaries implemented through the self-energy, and (2) the recursive Green's function method from the [RGF](https://github.com/chiaweihsu/RGF) repository; both are efficient for 1D systems and 2D systems with small-to-medium transverse widths. 

To compute the complete field profiles, simply omit the argument <code>C</code> or  <code>out</code>, or set it to <code>[]</code>.

Detailed usage of these functions are given in the documentation section (comment lines at the beginning) of the <code>.m</code> function files in the [src](./src) folder. One can, for example, type <code>help mesti</code> in MATLAB to see such documentation.

In addition, the user can use the function <code>mesti_build_channels()</code> to build the input and/or output matrices when using <code>mesti()</code>, or to determine which subset of the channels are of interest when using <code>mesti2s()</code>.

MESTI has two other general-purpose functions <code>mesti_build_fdfd_matrix()</code> and <code>mesti_matrix_solver()</code>; a typical user shouldn't need to directly use these two functions, but they are available.

## Examples

Several examples are given in the [examples](./examples) folder, which illustrate the usage and functionalities of MESTI. Each example has its own folder, with its <code>.m</code> script, auxiliary files specific to that example, and a <code>README.md</code> page that shows the code with its outputs:

1. [Fabry–Pérot etalon](./examples/1d_fabry_perot): 1D, using <code>mesti2s()</code>, with comparison to analytic solution.
2. [Distributed Bragg reflector](./examples/1d_distributed_bragg_reflector): 1D, using <code>mesti2s()</code>, with comparison to analytic solution.
3. [Open channel in a scattering medium](./examples/2d_open_channel_through_disorder): 2D, using <code>mesti2s()</code>, transmission matrix & field profile with customized wavefronts.
4.  [Reflection matrix in Gaussian-beam basis](./examples/2d_reflection_matrix_Gaussian_beams): 2D, using <code>mesti()</code>, reflection matrix in customized basis, for system open on all sides.
5. [Meta-atom design for metasurfaces](./examples/2d_meta_atom_design): 2D, using <code>mesti2s()</code> with Bloch periodic boundary for individual unit cells of a metasurface.
6. [Millimeter-wide metalens](./examples/2d_metalens_full): 2D, using <code>mesti()</code> with compressed input/output matrices (SCSA-c) for very-large-area simulations.

More examples will be added later.

## Gallery
Here are some images from the examples above:

1. Propagation through a Fabry–Pérot etalon
<img src="./examples/1d_fabry_perot/fabry_perot_field_profile.gif" width="336" height="264"> 

2. Open channel propagating through a strongly scattering medium
<img src="./examples/2d_open_channel_through_disorder/disorder_open_channel.gif" width="530" height="398"> 

3. Reflection matrix of a scatterer in Gaussian-beam basis:
<img src="./examples/2d_reflection_matrix_Gaussian_beams/reflection_matrix_Gaussian_beams.gif" width="447" height="257"> 

4. Full-wave simulation of a mm-wide hyperbolic metalens 
<img src="./examples/2d_metalens_full/metalens_intensity_profile_0_degree.png" width="400" height="300"> 

5. Angle dependence of the metalens above
<img src="./examples/2d_metalens_full/metalens_Strehl_ratio_and_transmission_eff.png" width="400" height="300"> 

## Reference

For more information on the theory, capability, and benchmarks (*e.g.*, scaling of computing time, memory usage, and accuracy), please see:

- Ho-Chun Lin, Zeyu Wang, and Chia Wei Hsu, "Full-wave solver for massively-multi-channel optics using Schur complement,"  [arXiv:2203.xxxxx](https://arxiv.org/abs/2203.xxxxx) (2022).

Please cite the above paper if you use MESTI for your work.
