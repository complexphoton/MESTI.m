function [S, stat] = mesti(syst, B, C, D, opts)
%MESTI Frequency-domain electromagnetics simulations.
%   [field_profiles, stat] = MESTI(syst, B) returns the spatial field profiles
%   E_z(x,y) for 2D transverse-magnetic (TM) waves satisfying
%      [- (d/dx)^2 - (d/dy)^2 - (omega/c)^2*epsilon_r(x,y)] E_z(x,y) = 0.
%   The frequency omega, the relative permittivity profile epsilon_r(x,y), and
%   the boundary conditions are aspecified by structure 'syst'. Each column of
%   matrix 'B' specifies a distinct input source profile. The returned
%   'field_profiles' is a 3D array, with field_profiles(:,:,i) being the field
%   profile given the i-th input source profile. The statistics of the
%   computation is returned in structure 'stat'.
%
%   MESTI uses finite-difference discretization, after which the differential
%   operator becomes an nx*ny-by-nx*ny sparse matrix A where [ny, nx] =
%   size(syst.epsilon), and the field profiles are given by
%   reshape(inv(A)*B, ny, nx, []).
%
%   [S, stat] = MESTI(syst, B, C) returns S = C*inv(A)*B where the solution
%   inv(A)*B is projected onto the output channels or locations of interest
%   throught matrix C; each row of matrix 'C' is a distinct output profile,
%   discretized into a 1-by-(nx*ny) vector in the same oreder as matrix A and
%   reshape(syst.epsilon, 1, []). When the MUMPS function zmumps() is
%   available, this is done by computing the Schur complement of an augmented
%   matrix K = [A,B;C,0] through a partial factorization.
%
%   [S, stat] = MESTI(syst, B, C, D) returns S = C*inv(A)*B - D. This can be
%   used for the computation of scattering matrices, where S is the scattering
%   matrix, and matrix D can be derived analytically or computed as D =
%   C*inv(A0)*B - S0 from a reference system A0 for which the scattering matrix
%   S0 is known.
%
%   [field_profiles, stat] = MESTI(syst, B, [], [], opts),
%   [S, stat] = MESTI(syst, B, C, [], opts), and
%   [S, stat] = MESTI(syst, B, C, D, opts) allow detailed options to be specified
%   with structure 'opts'.
%
%   This file checks and parses the parameters, and it can build matrices B and
%   C from the non-zero parts specified by the user (see details below). It
%   calls function mesti_build_fdfd_matrix() to build matrix A and function
%   mesti_matrix_solver() to compute C*inv(A)*B or inv(A)*B, where most of the
%   computation is done.
%
%   === Input Arguments ===
%   syst (scalar structure; required):
%      A structure that specifies the system, used to build the FDFD matrix A.
%      It contains the following fields:
%      syst.epsilon (numeric matrix, real or complex; required):
%         Discretized relative permittivity profile, with syst.epsilon(m,n)
%         being the relative permittivity epsilon_r(x,y) averaged over a square
%         with area (syst.dx)^2 centered at x = x_n = (n-0.5)*syst.dx, y = y_m =
%         (m-0.5)*syst.dx. The system size is size(syst.epsilon) = (ny, nx).
%         Note that y corresponds to the first index m, and x corresponds to the
%         second index n.
%            The positive imaginary part of syst.epsilon describes absorption,
%         and the negative imaginary part describes linear gain.
%      syst.length_unit (anything; optional):
%         Length unit, such as micron, nm, or some reference wavelength. This
%         code only uses dimensionless quantities, so syst.length_unit is never
%         used. This syst.length_unit is meant to help the user interpret the
%         units of (x,y), dx, wavelength, kx_B, ky_B, etc.
%      syst.wavelength (numeric scalar, real or complex; required):
%         Vacuum wavelength 2*pi*c/omega, in units of syst.length_unit.
%      syst.dx (positive scalar; required):
%         Discretization grid size, in units of syst.length_unit.
%      syst.PML (scalar structure or cell array; optional):
%         Parameters of the perfectly matched layer (PML) used to simulate an
%         open boundary. Note that PML is not a boundary condition; it is a
%         layer placed within the simulation domain (just before the boundary)
%         that attenuates outgoing waves with minimal reflection. In mesti(), the
%         PML starts from within syst.epsilon and ends at the first or last
%         pixel of syst.epsilon. (Note: this is different from the function
%         mesti2s() that handles two-sided geometries, where the homogeneous
%         spaces on the left and right are specified separately through
%         syst.epsilon_L and syst.epsilon_R, and where PML is placed in such
%         homogeneous space, outside of syst.epsilon.) When only one set of PML
%         parameters is used in the system (as is the most common), such
%         parameters can be specified by using a scalar structure syst.PML
%         containing the following fields:
%            npixels (positive integer scalar; required): Number of PML pixels.
%               Note this is within syst.epsilon.
%            direction (character vector; optional): Direction(s) where PML is
%               placed. Available choices are (case-insensitive):
%                  'all' - (default) PML in both x and y directions
%                  'x'   - PML in x direction
%                  'y'   - PML in y direction
%            side (character vector; optional): Side(s) where PML is placed.
%               Available choices are (case-insensitive):
%                  'all'  - (default) PML on both sides
%                  'low'  - one-sided PML; end at the first pixel (n=1 or m=1)
%                  'high' - one-sided PML; end at the last pixel (n=nx or m=ny)
%                  'both' - same as 'all'
%                  '-'    - same as 'low'
%                  '+'    - same as 'high'
%            power_sigma (non-negative scalar; optional): Power of the
%               polynomial grading for the conductivity sigma; defaults to 3.
%            sigma_max_over_omega (non-negative scalar; optional):
%               Conductivity at the end of the PML; defaults to
%                  0.8*(power_sigma+1)/((2*pi/wavelength)*dx*sqrt(epsilon_bg)).
%               where epsilon_bg is the average relative permittivity along the
%               last slice of syst.epsilon of the PML. This is used to attenuate
%               propagating waves.
%            power_kappa (non-negative scalar; optional): Power of the
%               polynomial grading for the real-coordinate-stretching factor
%               kappa; defaults to 3.
%            kappa_max (real scalar no smaller than 1; optional):
%               Real-coordinate-stretching factor at the end of the PML;
%               defaults to 15. This is used to accelerate the attenuation of
%               evanescent waves. kappa_max = 1 means no real-coordinate
%               stretching.
%            power_alpha (non-negative scalar; optional): Power of the
%               polynomial grading for the CFS alpha factor; defaults to 1.
%            alpha_max_over_omega (non-negative scalar; optional): Complex-
%               frequency-shifting (CFS) factor at the beginning of the PML.
%               This is typically used in time-domain simulations to suppress
%               late-time (low-frequency) reflections. We don't use it by
%               default (alpha_max_over_omega = 0) since we are in frequency
%               domain.
%         We use the following PML coordinate-stretching factor:
%            s(u) = kappa(u) + sigma(u)./(alpha(u) - i*omega)
%         with
%            sigma(u)/omega = sigma_max_over_omega*(u.^power_sigma),
%            kappa(u) = 1 + (kappa_max-1)*(u.^power_kappa),
%            alpha(u)/omega = alpha_max_over_omega*((1-u).^power_alpha),
%         where omega is frequency, and u goes linearly from 0 at the beginning
%         of the PML to 1 at the end of the PML. 
%            By default, syst.PML = {}, which means no PML on any side. PML is
%         only placed on the side(s) specified by syst.PML.
%            When multiple sets of PML parameters are used in the system (e.g.,
%         a thinner PML on one side, a thicker PML on another side), these
%         parameters can be specified with a cell array
%            syst.PML = {PML_1, PML_2, ...},
%         with PML_1 and PML_2 each being a structure containing the above
%         fields; they can specify different PML parameters on different sides.
%         Each side cannot be specified more than once.
%            With real-coordinate stretching, PML can attenuate evanescent waves
%         more efficiently than free space, so there is no need to place free
%         space in front of PML.
%            The PML thickness should be chosen based on the acceptable level of
%         reflectivity given the discretization resolution and the range of wave
%         numbers (i.e., angles) involved; more PML pixels gives lower
%         reflectivity. Typically 10-40 pixels are sufficient.
%      syst.PML_type (character vector; optional):
%         Type of PML. Available choices are (case-insensitive):
%            'UPML'   - (default) uniaxial PML
%            'SC-PML' - stretched-coordinate PML
%         The two are mathematically equivalent, but matrix A using UPML is
%         symmetric (unless Bloch periodic boundary is used) while that using
%         SC-PML has lower condition number.
%      syst.xBC (character vector; optional):
%         Boundary condition (BC) at the two ends in x direction, effectively
%         specifying E_z(m,n) at n=0 and n=nx+1 which are one pixel beyond our
%         computation domain. Available choices are (case-insensitive):
%            'Bloch'     - E_z(m,n+nx) = E_z(m,n)*exp(1i*syst.kx_B*nx*syst.dx)
%            'periodic'  - E_z(m,n+nx) = E_z(m,n)
%            'PEC'       - E_z(m,0) = E_z(m,nx+1) = 0
%            'PMC'       - E_z(m,0) = E_z(m,1); E_z(m,nx+1) = E_z(m,nx)
%            'PECPMC'    - E_z(m,0) = 0; E_z(m,nx+1) = E_z(m,nx)
%            'PMCPEC'    - E_z(m,0) = E_z(m,1); E_z(m,nx+1) = 0
%            'Dirichlet' - same as 'PEC'
%            'Neumann'   - same as 'PMC'
%         By default, syst.xBC = 'Bloch' if syst.kx_B is given; syst.xBC = 'PEC'
%         otherwise. The choice of syst.xBC has little effect on the numerical
%         accuracy when PML is used.
%      syst.kx_B (real scalar; optional):
%         Bloch wave number in x direction, in units of 1/syst.length_unit.
%         syst.kx_B is only used when syst.xBC = 'Bloch'. It is allowed to
%         specify a complex-valued syst.kx_B, but a warning will be displayed.
%      syst.yBC (character vector; optional):
%         Boundary condition in y direction, analogous to syst.xBC.
%      syst.ky_B (real scalar; optional):
%         Bloch wave number in y direction, analogous to syst.kx_B.
%      syst.self_energy (sparse matrix; optional):
%         Self-energy matrix, used as A = A - syst.self_energy to achieve exact
%         radiation boundary condition. In mesti2s() when syst.xBC = 'outgoing',
%         the self-energy matrix will be built and passed to mesti(). On the
%         sides where self-energy is used, Dirichlet boundary condition should
%         be use with no PML.
%   B (numeric matrix or structure array; required):
%      Matrix specifying the input in the C*inv(A)*B - D or C*inv(A)*B or
%      inv(A)*B returned. When the input argument B is a matrix, it is directly
%      used, and size(B,1) must equal ny*nx = numel(syst.epsilon); different
%      columns of B correspond to different inputs.
%         Instead of specifying matrix B directly, one can specify only its
%      non-zero parts, from which mesti() will build the sparse matrix B. To do
%      so, B in the input argument should be set as a structure array; here we
%      refer to such structure array as B_struct to distinguish it from the
%      resulting matrix B. If for every column of matrix B, all of its non-zero
%      elements are spatially located within a rectangle (e.g., line sources or
%      block sources), one can use the following fields:
%         B_struct.pos (four-element integer vector): B_struct.pos =
%            [m1, n1, h, w] specifies the location and the size of the
%            rectangle. Here, (m1, n1) is the index of the (y,x) coordinate of
%            the smaller-y, smaller-x corner of the rectangle, at the location
%            of syst.epsilon(m1, n1); (h, w) is the height and width of the
%            rectangle, such that (m2, n2) = (m1+h-1, n1+w-1) is the index of
%            the higher-index corner of the rectangle, at the location of
%            syst.epsilon(m2, n2).
%         B_struct.data (2D or 3D numeric array): non-zero elements of matrix B
%            within the rectangle specified by B_struct.pos. When it is a 3D
%            array, B_struct.data(m',n',a) is the a-th input source at the
%            location of syst.epsilon(m=m1+m'-1, n=n1+n'-1), which becomes
%            B(m+(n-1)*ny, a). In other words, B_struct.data(:,:,a) gives the
%            sources at the rectangle syst.epsilon(m1+(0:(h-1)), n1+(0:(w-1))).
%            So, size(B_struct.data, [1,2]) must equal [h, w], and
%            size(B_struct.data, 3) is the number of inputs. Alternatively,
%            B_struct.data can be a 2D array that is equivalent to
%            reshape(data_in_3D_array, h*w, []), in which case
%            size(B_struct.data, 2) is the number of inputs; in this case,
%            B_struct.data can be a sparse matrix, and its sparsity will be
%            preserved when building matrix B.
%         If different inputs are located within different rectangles (e.g.,
%      inputs from line sources on the left and separate inputs from line
%      sources on the right), B_struct can be a structure array with multiple
%      elements [e.g., B_struct(1).pos and B_struct(1).data specify line sources
%      on the left; B_struct(2).pos and B_struct(2).data specify line sources on
%      the left]; these inputs are treated separately, and the total number of
%      inputs is size(B_struct(1).data, 3) + size(B_struct(2).data, 3) + ... +
%      size(B_struct(end).data, 3).
%         If the non-zero elements of matrix B do not have rectangular shapes in
%      space [e.g., for total-field/scattered-field (TF/SF) simulations], one
%      can use a structure array with the following fields:
%         B_struct.ind (integer vector): linear indices of the spatial
%            locations of the non-zero elements of matrix B, such that
%            syst.epsilon(B_struct.ind) are the points where the source is
%            placed. Such linear indices can be constructed from sub2ind().
%         B_struct.data (2D numeric matrix): non-zero elements of matrix B at
%            the locations specified by B_struct.ind. Specifically,
%            B_struct.data(i,a) is the a-th input source at the location of
%            syst.epsilon(B_struct.ind(i)), which becomes B(B_struct.ind(i), a).
%            So, size(B_struct.data, 1) must equal numel(B_struct.ind), and
%            size(B_struct.data, 2) is the number of inputs.
%         Similarly, one can use B_struct(1).ind, B_struct(2).ind etc together
%      with B_struct(1).data, B_struct(2).data etc to specify inputs at
%      different sets of locations. Every element of the structure array must
%      have the same fields [e.g., one cannot specify B_struct(1).pos and
%      B_struct(2).ind], so the more general B_struct.ind syntax should be used
%      when some of the inputs are rectangular and some are not.
%   C (numeric matrix or structure array or 'transpose(B)' or []; optional):
%      Matrix specifying the output projections in the C*inv(A)*B - D or
%      C*inv(A)*B returned. When the input argument C is a matrix, it is
%      directly used, and size(C,2) must equal ny*nx = numel(syst.epsilon);
%      different rows of C correspond to different outputs.
%         Scattering matrix computations often have C = transpose(B); if that
%      is the case, the user can set C = 'transpose(B)' as a character vector,
%      and it will be replaced by transpose(B) in the code. Doing so has an
%      advantage: if matrix A is symmetric (which is the case with UPML without
%      Bloch periodic boundary), C = 'transpose(B)', opts.solver = 'MUMPS', and
%      opts.method = 'APF', the matrix K = [A,B;C,0] will be treated as
%      symmetric when computing its Schur complement to lower computing time and
%      memory usage.
%         For field profile computations, the user can simply omit C from the
%      input argument, as in mesti(syst, B), if there is no need to change the
%      default opts. If opts is needed, the user can use
%      mesti(syst, B, [], [], opts), namely setting C = [] and D = [].
%         Similar to B, here one can specify only the non-zero parts of the
%      output matrix C, from which mesti() will build the sparse matrix C. The
%      syntax is the same as for B, summarized below. If for every row of matrix
%      C, all of its non-zero elements are spatially located withing a rectangle
%      (e.g., projection of fields on a line), one can set the input argument C
%      to be a structure array (referred to as C_struct below) with the
%      following fields:
%         C_struct.pos (four-element integer vector): C_struct.pos =
%            [m1, n1, h, w] specifies the location and the size of the
%            rectangle. Here, (m1, n1) is the index of the (y,x) coordinate of
%            the smaller-y, smaller-x corner of the rectangle, at the location
%            of syst.epsilon(m1, n1); (h, w) is the height and width of the
%            rectangle, such that (m2, n2) = (m1+h-1, n1+w-1) is the index of
%            the higher-index corner of the rectangle, at the location of
%            syst.epsilon(m2, n2).
%         C_struct.data (2D or 3D numeric array): non-zero elements of matrix C
%            within the rectangle specified by C_struct.pos. When it is a 3D
%            array, C_struct.data(m',n',b) is the b-th output projection at the
%            location of syst.epsilon(m=m1+m'-1, n=n1+n'-1), which becomes
%            C(b, m+(n-1)*ny). In other words, C_struct.data(:,:,b) gives the
%            projection at the rectangle syst.epsilon(m1+(0:(h-1)), n1+(0:(w-1))).
%            So, size(C_struct.data, [1,2]) must equal [h, w], and
%            size(C_struct.data, 3) is the number of outputs. Alternatively,
%            C_struct.data can be a 2D array that is equivalent to
%            reshape(data_in_3D_array, h*w, []), in which case
%            size(C_struct.data, 2) is the number of outputs; in this case,
%            C_struct.data can be a sparse matrix, and its sparsity will be
%            preserved when building matrix C.
%         If the non-zero elements of matrix C do not have rectangular shapes in
%      space [e.g., for near-field-to-far-field transformations], one can set C
%      to a structure array with the following fields:
%         C_struct.ind (integer vector): linear indices of the spatial
%            locations of the non-zero elements of matrix C, such that
%            syst.epsilon(C_struct.ind) are the points where the projection is
%            placed. Such linear indices can be constructed from sub2ind().
%         C_struct.data (2D numeric matrix): non-zero elements of matrix C at
%            the locations specified by C_struct.ind. Specifically,
%            C_struct.data(i,b) is the b-th projection at the location of
%            syst.epsilon(C_struct.ind(i)), which becomes C(b, C_struct.ind(i)).
%            So, size(C_struct.data, 1) must equal numel(C_struct.ind), and
%            size(C_struct.data, 2) is the number of outputs.
%         Like in B_struct, one can use structure arrays with multiple elements
%      to specify outputs at different spatial locations.
%   D (numeric matrix or []; optional):
%      Matrix specifying the background in the absence of scatterers for
%      scattering matrix computations in the form of C*inv(A)*B - D; size(D,1)
%      must equal size(C,1), and size(D,2) must equal size(B,2).
%         When D = [], it will not be subtracted from C*inv(A)*B. For field
%      profile computations where C = [], the user must also set D = [].
%   opts (scalar structure; optional, defaults to an empty struct):
%      A structure that specifies the options of computation; defaults to an
%      empty structure. It can contain the following fields (all optional):
%      opts.verbal (logical scalar; optional, defaults to true):
%         Whether to print info and timing to the standard output.
%      opts.prefactor (numeric scalar, real or complex; optional):
%         When opts.prefactor is given, mesti() will return
%         opts.prefactor*C*inv(A)*B - D or opts.prefactor*C*inv(A)*B or
%         opts.prefactor*inv(A)*B. Such prefactor makes it easier to use C =
%         transpose(B) to take advantage of reciprocity. Defaults to 1.
%      opts.solver (character vector; optional):
%         The software used for sparse matrix factorization. Available choices
%         are (case-insensitive):
%            'MUMPS'  - (default) Uses MUMPS. Its MATLAB interface zmumps.m must
%                       be in MATLAB's search path. This is much faster and uses
%                       less memory.
%            'MATLAB' - Uses the built-in lu() function in MATLAB, which uses
%                       UMFPACK with AMD ordering. This requires no installation
%                       but is much slower. This is be used by default if
%                       zmumps.m is not found in the search path.
%      opts.method (character vector; optional):
%         The solution method. Available choices are (case-insensitive):
%            'APF' - Augmented partial factorization. When opts.solver =
%                    'MUMPS', C*inv(A)*B is obtained through the Schur
%                    complement of an augmented matrix K = [A,B;C,0] using a
%                    partial factorization; this is the true APF. When
%                    opts.solver = 'MATLAB', C*inv(A)*B is obtained as
%                    C*inv(U)*inv(L)*B with optimized grouping, which is not the
%                    true APF but is slightly better than factorize_and_solve.
%                    Cannot be used for computing the full field profile
%                    inv(A)*B or with iterative refinement.
%            'FS'  - Factorize and solve. Factorize A=L*U, solve for inv(A)*B
%                    with forward and backward substitutions, and optionally
%                    project with C.
%            'factorize_and_solve' - Same as 'FS'.
%         By default, if C is given and opts.iterative_refinement = false, then
%         'APF' is used. Otherwise, 'factorize_and_solve' is used.
%      opts.clear_BC (logical scalar; optional, defaults to false):
%         When opts.clear_BC = true, variables 'B' and 'C' will be cleared in
%         the caller's workspace to reduce peak memory usage. Can be used when B
%         and/or C take up significant memory and are not needed after calling
%         mesti().
%      opts.clear_syst (logical scalar; optional, defaults to false):
%         When opts.clear_syst = true, variable 'syst' will be cleared in the
%         caller's workspace to reduce peak memory usage. Can be used when
%         syst.epsilon and/or syst.self_energy take up significant memory and
%         are not needed after calling mesti().
%      opts.clear_memory (logical scalar; optional, defaults to true):
%         Whether or not to clear variables inside mesti() to reduce peak memory
%         usage.
%      opts.verbal_solver (logical scalar; optional, defaults to false):
%         Whether to have the solver print detailed information to the standard
%         output. Note the behavior of output from MUMPS depends on compiler.
%      opts.use_METIS (logical scalar; optional, defaults to false):
%         Whether to use METIS (instead of the default AMD) to compute the
%         ordering in MUMPS. Using METIS can sometimes reduce memory usage
%         and/or factorization and solve time, but it typically takes longer at
%         the analysis (i.e., ordering) stage.
%      opts.nrhs (positive integer scalar; optional):
%         The number of right-hand sides (number of columns of the input matrix
%         B) to consider simultaneously, used only when opts.method =
%         'factorize_and_solve' and C is given. Defaults to 1 if
%         opts.iterative_refinement = true, 10 if opts.solver = 'MUMPS' with
%         opts.iterative_refinement = false, 4 otherwise.
%      opts.store_ordering (logical scalar; optional, defaults to false):
%         Whether to store the ordering sequence (permutation) for matrix A or
%         matrix K; only possible when opts.solver = 'MUMPS'. If
%         opts.store_ordering = true, the ordering will be returned in
%         stat.ordering.
%      opts.ordering (positive integer vector; optional):
%         A user-specified ordering sequence for matrix A or matrix K, used only
%         when opts.solver = 'MUMPS'. Using the ordering from a previous
%         computation can speed up the analysis stage, but the matrix size must
%         be the same.
%      opts.nthreads_OMP (positive integer scalar; optional):
%         Number of OpenMP threads used in MUMPS; overwrites the OMP_NUM_THREADS
%         environment variable.
%      opts.iterative_refinement (logical scalar; optional, defaults to false):
%         Whether to use iterative refinement in MUMPS to lower round-off
%         errors. Iterative refinement can only be used when opts.solver =
%         'MUMPS' and opts.method = 'factorize_and_solve' and C is given, in
%         case opts.nrhs must equal 1. When iterative refinement is used, the
%         relevant information will be returned in stat.itr_ref_nsteps,
%         stat.itr_ref_omega_1, and stat.itr_ref_omega_2.
%
%   === Output Arguments ===
%   S (full numeric matrix or 3d array):
%      C*inv(A)*B or reshape(inv(A)*B, ny, nx, []).
%   stat (scalar structure):
%      A structure that contains the following fields:
%      stat.opts (scalar structure):
%         The final 'opts' used, excluding the user-specified matrix ordering.
%      stat.timing (scalar structure):
%         A structure containing timing of the various stages, in seconds, in
%         fields 'total', 'init', 'build', 'analyze', 'factorize', 'solve'.
%      stat.xPML (two-element cell array; optional);
%         PML parameters on the low and high sides of x direction, if used.
%      stat.yPML (two-element cell array; optional);
%         PML parameters on the low and high sides of y direction, if used.
%      stat.ordering_method (character vector; optional):
%         Ordering method used in MUMPS.
%      stat.ordering (positive integer vector; optional):
%         Ordering sequence returned by MUMPS when opts.store_ordering = true.
%      stat.itr_ref_nsteps (integer vector; optional):
%         Number of steps of iterative refinement for each input, if
%         opts.iterative_refinement = true; 0 means no iterative refinement.
%      stat.itr_ref_omega_1 (real vector; optional):
%         Scaled residual omega_1 at the end of iterative refinement for each
%         input; see MUMPS user guide section 3.3.2 for definition.
%      stat.itr_ref_omega_2 (real vector; optional):
%         Scaled residual omega_2 at the end of iterative refinement for each
%         input; see MUMPS user guide section 3.3.2 for definition.
%
%   See also: mesti_build_fdfd_matrix, mesti_matrix_solver, mesti2s

%% Part 1.1: Check validity of syst, assign default values to its fields, and parse BC and PML specifications

t0 = clock;

if nargin < 2; error('Not enough input arguments.'); end
if ~(isstruct(syst) && isscalar(syst)); error('Input argument ''syst'' must be a scalar structure.'); end
if ~isfield(syst, 'epsilon');           error('Input argument ''syst'' must have field ''epsilon''.'); end
if ~isfield(syst, 'wavelength');        error('Input argument ''syst'' must have field ''wavelength''.'); end
if ~isfield(syst, 'dx');                error('Input argument ''syst'' must have field ''dx''.'); end
if ~(isnumeric(syst.epsilon)    && ismatrix(syst.epsilon));    error('syst.epsilon must be a numeric matrix.'); end
if ~(isnumeric(syst.wavelength) && isscalar(syst.wavelength)); error('syst.wavelength must be a numeric scalar.'); end
if ~(isreal(syst.dx) && isscalar(syst.dx) && syst.dx > 0);     error('syst.dx must be a positive scalar.'); end

% Check that the user did not accidentally use options only in mesti2s()
if isfield(syst, 'epsilon_L') && ~isempty(syst.epsilon_L)
    warning('syst.epsilon_L is not used in mesti(); will be ignored.');
end
if isfield(syst, 'epsilon_R') && ~isempty(syst.epsilon_R)
    warning('syst.epsilon_R is not used in mesti(); will be ignored.');
end

% Number of grid points in y and x
[ny, nx] = size(syst.epsilon);
nxy = nx*ny;

% Check boundary condition in x
if isfield(syst, 'kx_B') && ~isempty(syst.kx_B)
    if ~(isnumeric(syst.kx_B) && isscalar(syst.kx_B))
        error('syst.kx_B must be a numeric scalar, if given.');
    elseif (isfield(syst, 'xBC') && ~isempty(syst.xBC)) && (iscell(syst.xBC) || ~strcmpi(syst.xBC, 'Bloch'))
        error('When syst.kx_B is given, syst.xBC must be ''Bloch'' if specified.');
    end
    syst.xBC = 'Bloch';
    % mesti_build_fdfd_matrix() uses (kx_B,ky_B)*periodicity as the input arguments xBC and yBC for Bloch BC
    xBC = (syst.kx_B)*(nx*syst.dx); % dimensionless
else
    % Defaults to PEC boundary condition unless syst.kx_B is given
    if ~isfield(syst, 'xBC') || isempty(syst.xBC)
        syst.xBC = 'PEC';
    elseif ~((ischar(syst.xBC) && isrow(syst.xBC)) || (isstring(syst.xBC) && isscalar(syst.xBC)))
        error('syst.xBC must be a character vector or string, if given.');
    elseif ~ismember(lower(syst.xBC), {'bloch', 'periodic', 'dirichlet', 'pec', 'neumann', 'pmc', 'dirichletneumann', 'pecpmc', 'neumanndirichlet', 'pmcpec'})
        error('syst.xBC = ''%s'' is not a supported option; type ''help mesti'' for supported options.', syst.xBC);
    elseif strcmpi(syst.xBC, 'Bloch')
        error('syst.xBC = ''Bloch'' but syst.kx_B is not given.');
    end
    xBC = syst.xBC;
end
if isfield(syst, 'kx') && ~isempty(syst.kx)
    warning('syst.kx will not be used; use syst.kx_B for Bloch wave number in x.');
end

% Check boundary condition in y
if isfield(syst, 'ky_B') && ~isempty(syst.ky_B)
    if ~(isnumeric(syst.ky_B) && isscalar(syst.ky_B))
        error('syst.ky_B must be a numeric scalar, if given.');
    elseif (isfield(syst, 'yBC') && ~isempty(syst.yBC)) && (iscell(syst.yBC) || ~strcmpi(syst.yBC, 'Bloch'))
        error('When syst.ky_B is given, syst.yBC must be ''Bloch'' if specified.');
    end
    syst.yBC = 'Bloch';
    % mesti_build_fdfd_matrix() uses (kx_B,ky_B)*periodicity as the input arguments xBC and yBC for Bloch BC
    yBC = (syst.ky_B)*(ny*syst.dx); % dimensionless
else
    % Defaults to PEC boundary condition unless syst.ky_B is given
    if ~isfield(syst, 'yBC') || isempty(syst.yBC)
        syst.yBC = 'PEC';
    elseif ~((ischar(syst.yBC) && isrow(syst.yBC)) || (isstring(syst.yBC) && isscalar(syst.yBC)))
        error('syst.yBC must be a character vector or string, if given.');
    elseif ~ismember(lower(syst.yBC), {'bloch', 'periodic', 'dirichlet', 'pec', 'neumann', 'pmc', 'dirichletneumann', 'pecpmc', 'neumanndirichlet', 'pmcpec'})
        error('syst.yBC = ''%s'' is not a supported option; type ''help mesti'' for supported options.', syst.yBC);
    elseif strcmpi(syst.yBC, 'Bloch')
        error('syst.yBC = ''Bloch'' but syst.ky_B is not given.');
    end
    yBC = syst.yBC;
end
if isfield(syst, 'ky') && ~isempty(syst.ky)
    warning('syst.ky will not be used; use syst.ky_B for Bloch wave number in y.');
end

% Defaults to no PML anywhere
if ~isfield(syst, 'PML') || isempty(syst.PML)
    syst.PML = {};
elseif ~((isstruct(syst.PML) && isscalar(syst.PML)) || iscell(syst.PML))
    error('syst.PML must be a scalar structure or a cell array, if given.');
elseif isstruct(syst.PML)
    % convert to a single-element cell array if only one set of PML spec is given
    syst.PML = {syst.PML};
end

% Parse the user-specified PML parameters to PML on the four sides
% PML_list = {xPML_low, xPML_high, yPML_low, yPML_high}
PML_list = {[], [], [], []};
str_sides = {'-x', '+x', '-y', '+y'};
use_PML = false;
for ii = 1:numel(syst.PML)
    use_PML = true;
    PML_ii = syst.PML{ii};
    if ~(isstruct(PML_ii) && isscalar(PML_ii))
        error('syst.PML{%d} must be a scalar structure.', ii)
    end

    % Number of PML pixels must be given
    % Other fields are optional and will be checked in mesti_build_fdfd_matrix()
    if ~isfield(PML_ii, 'npixels') || isempty(PML_ii.npixels)
        error('syst.PML{%d} must contain field ''npixels''.', ii);
    end

    % If PML is specified, we put it on both x and y directions by default
    if ~isfield(PML_ii, 'direction') || isempty(PML_ii.direction)
        PML_ii.direction = 'all';
    elseif ~((ischar(PML_ii.direction) && isrow(PML_ii.direction)) || (isstring(PML_ii.direction) && isscalar(PML_ii.direction)))
        error('syst.PML{%d}.direction must be a character vector or string, if given.', ii);
    elseif ~ismember(lower(PML_ii.direction), {'all', 'x', 'y'})
        error('syst.PML{%d}.direction = ''%s'' is not a supported option; use ''all'', ''x'', or ''y''.', ii, PML_ii.direction);
    end

    % If PML is specified, we put it on both low and high sides by default
    if ~isfield(PML_ii, 'side') || isempty(PML_ii.side)
        PML_ii.side = 'all';
    elseif ~((ischar(PML_ii.side) && isrow(PML_ii.side)) || (isstring(PML_ii.side) && isscalar(PML_ii.side)))
        error('syst.PML{%d}.side must be a character vector or string, if given.', ii);
    elseif ~ismember(lower(PML_ii.side), {'all', 'both', 'low', 'high', '-', '+'})
        error('syst.PML{%d}.side = ''%s'' is not a supported option; use ''all'', ''low'', or ''high''.', ii, PML_ii.side);
    end

    % Convert {PML_ii.direction and PML_ii.side} to a list of the PML locations
    % 1=xPML_low, 2=xPML_high, 3=yPML_low, 4=yPML_high
    if strcmpi(PML_ii.direction, 'all') % x & y
        if strcmpi(PML_ii.side, 'all') || strcmpi(PML_ii.side, 'both') % low & high
            ind_ii = [1,2,3,4];
        elseif strcmpi(PML_ii.side, 'low') || strcmpi(PML_ii.side, '-')
            ind_ii = [1,3];
        else % PML_ii.side = 'high' or '+'
            ind_ii = [2,4];
        end
    elseif strcmpi(PML_ii.direction, 'x')
        if strcmpi(PML_ii.side, 'all') || strcmpi(PML_ii.side, 'both') % low & high
            ind_ii = [1,2];
        elseif strcmpi(PML_ii.side, 'low') || strcmpi(PML_ii.side, '-')
            ind_ii = 1;
        else % PML_ii.side = 'high' or '+'
            ind_ii = 2;
        end
    else % PML_ii.direction = 'y'
        if strcmpi(PML_ii.side, 'all') || strcmpi(PML_ii.side, 'both') % low & high
            ind_ii = [3,4];
        elseif strcmpi(PML_ii.side, 'low') || strcmpi(PML_ii.side, '-')
            ind_ii = 3;
        else % PML_ii.side = 'high' or '+'
            ind_ii = 4;
        end
    end
    % These two fields are no longer needed
    PML_ii = rmfield(PML_ii, {'direction', 'side'});

    % Specify PML at those locations
    for jj = 1:numel(ind_ii)
        ind_side = ind_ii(jj);
        % Check that PML has not been specified at that location yet
        if ~isempty(PML_list{ind_side})
            error('PML on %s side is specified more than once in syst.PML.', str_sides{ind_side});
        end
        PML_list{ind_side} = PML_ii;
    end
end

% Convert to two separate cell arrays for mesti_build_fdfd_matrix()
xPML = PML_list(1:2); % {xPML_low, xPML_high}
yPML = PML_list(3:4); % {yPML_low, yPML_high}

% Use UPML by default as it produces a symmetric matrix A (unless Bloch periodic boundary is used)
if ~isfield(syst, 'PML_type') || isempty(syst.PML_type)
    syst.PML_type = 'UPML';
elseif ~((ischar(syst.PML_type) && isrow(syst.PML_type)) || (isstring(syst.PML_type) && isscalar(syst.PML_type)))
    error('syst.PML_type must be a character vector or string, if given.');
elseif ~ismember(lower(syst.PML_type), {'upml', 'sc-pml', 'scpml'})
    error('syst.PML_type = ''%s'' is not a supported option; use ''UPML'' or ''SC-PML''.', syst.PML_type);
end
if strcmpi(syst.PML_type, 'UPML')
    use_UPML = true;
else
    use_UPML = false;
end

%% Part 1.2: Check validity of the other input arguments and assign default values

if ~((ismatrix(B) && isnumeric(B)) || (isstruct(B) && ~isempty(B)))
    error('Input argument ''B'' must be a numeric matrix or a non-empty structure array.'); 
end

% C is an optional argument
if nargin < 3
    C = [];
end

% D is an optional argument
if nargin < 4
    D = [];
end
if ~((ismatrix(D) && isnumeric(D)) || isempty(D))
    error('Input argument ''D'' must be a numeric matrix or [], if given.');
end

% opts is an optional argument
if nargin < 5 || isempty(opts)
    opts = struct();
end
if ~(isstruct(opts) && isscalar(opts))
    error('Input argument ''opts'' must be a scalar structure or [], if given.');
end

if isempty(C) && ~isstruct(C)
    opts.return_field_profile = true;
elseif (ismatrix(C) && isnumeric(C)) || (isstruct(C) && ~isempty(C))
    opts.return_field_profile = false;
    use_transpose_B = false;
elseif isequal(C, 'transpose(B)')
    opts.return_field_profile = false;
    use_transpose_B = true;
else
    error('Input argument ''C'' must be a numeric matrix or a non-empty structure array or ''transpose(B)'' or [], if given.');
end

% Check that the user did not accidentally use options only in mesti2s()
if isfield(opts, 'symmetrize_K') && ~isempty(opts.symmetrize_K)
    error('opts.symmetrize_K is not used in mesti(); to symmetrize matrix K = [A,B;C,0], set C = ''transpose(B)'', make sure matrix A is symmetric (syst.PML_type = ''UPML'' and no Bloch periodic boundary), set opts.solver = ''MUMPS'', and set opts.method = ''APF''.');
end

% Turn on verbal output by default
if ~isfield(opts, 'verbal') || isempty(opts.verbal)
    opts.verbal = true;
elseif ~(islogical(opts.verbal) && isscalar(opts.verbal))
    error('opts.verbal must be a logical scalar, if given.');
end

% Defaults the prefactor to 1
if ~isfield(opts, 'prefactor') || isempty(opts.prefactor)
    opts.prefactor = 1;
elseif ~(isnumeric(opts.prefactor) && isscalar(opts.prefactor))
    error('opts.prefactor must be a numeric scalar, if given.');
end

% By default, we don't clear syst in the caller's workspace
if ~isfield(opts, 'clear_syst') || isempty(opts.clear_syst)
    opts.clear_syst = false;
elseif ~(islogical(opts.clear_syst) && isscalar(opts.clear_syst))
    error('opts.clear_syst must be a logical scalar, if given.');
end

% By default, we don't clear B and C in the caller's workspace
if ~isfield(opts, 'clear_BC') || isempty(opts.clear_BC)
    opts.clear_BC = false;
elseif ~(islogical(opts.clear_BC) && isscalar(opts.clear_BC))
    error('opts.clear_BC must be a logical scalar, if given.');
end

% By default, we will clear internal variables; only useed in mesti_matrix_solver()
% Note that mesti_matrix_solver() defaults opts.clear_memory to false because some users that use mesti_matrix_solver() directly may want to keep the input arguments A,B,C after calling it. But the opts.clear_memory in mesti() here only deals with the variables internal to mesti() and mesti_matrix_solver(); it doesn't deal with the input arguments provided by the user (which are specified by opts.clear_syst and opts.clear_BC), so it is safe to default opts.clear_memory to true here.
if ~isfield(opts, 'clear_memory') || isempty(opts.clear_memory)
    opts.clear_memory = true;
elseif ~(islogical(opts.clear_memory) && isscalar(opts.clear_memory))
    error('opts.clear_memory must be a logical scalar, if given.');
end

% The following fields of opts are not used in mesti() and will be checked/initialized in mesti_matrix_solver():
%    opts.solver
%    opts.method
%    opts.verbal_solver
%    opts.use_METIS
%    opts.nrhs
%    opts.store_ordering
%    opts.ordering
%    opts.nthreads_OMP
%    opts.iterative_refinement

if opts.verbal
    % print basic system info if the calling function is not mesti2s()
    st = dbstack;
    if numel(st) > 1 && strcmp(st(2).name,'mesti2s')
        called_from_mesti2s = true;
        fprintf('            ... ');
    else
        called_from_mesti2s = false;
        fprintf('System size: ny = %d, nx = %d\n', ny, nx);
        if use_PML
            fprintf('%s on ', syst.PML_type);
            for ind_side = 1:4
                if ~isempty(PML_list{ind_side})
                    fprintf('%s ', str_sides{ind_side});
                end
            end
            fprintf('sides; ');
        else
            fprintf('no PML; ');
        end
        fprintf('xBC = %s', syst.xBC);
        if strcmpi(syst.xBC, 'Bloch'); fprintf(' (kx_B = %.4f)', syst.kx_B); end
        fprintf('; yBC = %s', syst.yBC);
        if strcmpi(syst.yBC, 'Bloch'); fprintf(' (ky_B = %.4f)', syst.ky_B); end
        if isfield(syst, 'self_energy') && ~isempty(syst.self_energy); fprintf('; with self-energy'); end
        fprintf('\nBuilding B,C... ');
    end
end

t1 = clock; timing_init = etime(t1,t0); % Initialization time

%% Part 2.1: Build matrices B and C

% Build the input matrix B from its non-zero elements specified by user
if isstruct(B)
    B_struct = B;
    if ~isfield(B_struct, 'data')
        error('Input argument ''B'' must have field ''data'' when B is a structure array.');
    elseif ~isfield(B_struct, 'pos') && ~isfield(B_struct, 'ind')
        error('Input argument ''B'' must have field ''pos'' or ''ind'' when B is a structure array.');
    elseif isfield(B_struct, 'pos') && isfield(B_struct, 'ind')
        error('Input argument ''B'' cannot have both field ''pos'' and field ''ind'' when B is a structure array.');
    end
    if isfield(B_struct, 'pos')
        for ii = 1:numel(B_struct)
            pos = B_struct(ii).pos;
            if ~(isreal(pos) && isvector(pos) && numel(pos)==4 && isequal(pos, round(pos)) && min(pos)>0)
                error('B(%d).pos must be a positive integer vector with 4 elements.', ii);
            end
        end
    end

    % We first pick the most efficient way to build matrix B.
    % If all of the following are satisfied: (1) numel(B_struct) is small, (2) B_struct.pos is used, and (3) the rectangle specified by B_struct.pos is a single vertical slice or if it spans the full height of syst.epsilon, then we will stack reshaped B_struct(ii).data with zeros. This avoids the overhead of building B with index-value pairs.
    % If any of the above is not satisfied, we will build B with index-value pairs.
    use_iv_pairs = false;
    if numel(B_struct) > 10
        use_iv_pairs = true;
    elseif ~isfield(B_struct, 'pos')
        use_iv_pairs = true;
    else
        for ii = 1:numel(B_struct)
            pos = B_struct(ii).pos;
            if ~(pos(4) == 1 || pos(3) == ny)
                use_iv_pairs = true;
            end
        end
    end

    if use_iv_pairs
        % Construct matrix B from the complete set of index-value pairs
        N_tot = 0; % total number of non-zero elements in B
        for ii = 1:numel(B_struct)
            N_tot = N_tot + numel(B_struct(ii).data);
        end
        ind_list = zeros(N_tot, 0);
        a_list = zeros(N_tot, 0);
        v_list = zeros(N_tot, 0);
        N = 0;
        M = 0;
    else
        % Build matrix B incrementally
        B = sparse(nxy, 0);
    end

    % Loop over the structure array
    for ii = 1:numel(B_struct)
        data = B_struct(ii).data;
        if isfield(B_struct, 'pos')
            % B_struct(ii).pos specifies a rectangle inside syst.epsilon; (m1,n1) and (m2,n2) are its two diagonal corners
            pos = B_struct(ii).pos;
            m1 = pos(1); % first index in y
            n1 = pos(2); % first index in x
            m2 = m1 + pos(3) - 1; % last index in y
            n2 = n1 + pos(4) - 1; % last index in x
            nxy_data = pos(3)*pos(4); % number of elements in this rectangle
            if m1 > ny
                error('B(%d).pos(1) = %d exceeds ny = size(syst.epsilon, 1) = %d.', ii, m1, ny);
            elseif n1 > nx
                error('B(%d).pos(2) = %d exceeds nx = size(syst.epsilon, 2) = %d.', ii, n1, nx);
            elseif m2 > ny
                error('B(%d).pos(1) + B(%d).pos(3) - 1 = %d exceeds ny = size(syst.epsilon, 1) = %d.', ii, ii, m2, ny);
            elseif n2 > nx
                error('B(%d).pos(2) + B(%d).pos(4) - 1 = %d exceeds nx = size(syst.epsilon, 2) = %d.', ii, ii, n2, nx);
            elseif ~(ndims(data) <= 3 && isnumeric(data))
                error('B(%d).data must be a 2D or 3D numeric array when B.pos is given.', ii);
            elseif isequal(size(data, [1,2]), [pos(3), pos(4)]) % this way to use size() is supported starting MATLAB R2019b
                M_ii = size(data, 3); % number of inputs
            elseif size(data, 1) == nxy_data
                M_ii = size(data, 2); % number of inputs
            else
                error('size(B(%d).data) = [%d, %d, %d] is not compatible with [B(%d).pos(3), B(%d).pos(4)] = [%d, %d].', ii, size(data, 1), size(data, 2), size(data, 3), ii, ii, pos(3), pos(4));
            end
            if use_iv_pairs
                % convert to linear indices
                m_list = repmat((m1:m2).', 1, pos(4));
                n_list = repmat( n1:n2   , pos(3), 1);
                ind = reshape(sub2ind([ny, nx], m_list, n_list), nxy_data, 1);
            end
        else
            % B_struct(ii).ind specifies linear indices of syst.epsilon
            ind = B_struct(ii).ind;
            if ~(isreal(ind) && isvector(ind) && isequal(ind, round(ind)) && min(ind)>0 && max(ind)<=nxy && numel(ind)==numel(unique(ind)))
                error('B(%d).ind must be a vector with no repeated elements, where every element is an integer between 1 and nx*ny = %d.', ii, nxy);
            elseif ~(ismatrix(data) && isnumeric(data))
                error('B(%d).data must be a 2D numeric matrix when B.ind is given.', ii);
            elseif size(data, 1) ~= numel(ind)
                error('size(B(%d).data, 1) = %d does not match numel(B(%d).ind) = %d', ii, size(data, 1), ii, numel(ind));
            elseif ~iscolumn(ind)
                ind = ind(:);
            end
            nxy_data = numel(ind); % number of linear indices
            M_ii = size(data, 2); % number of inputs
        end

        if use_iv_pairs
            % Build index-value pairs: (ind_list, a_list, v_list)
            N_ii = nxy_data*M_ii; % number of non-zero elements in the ii-th part of matrix B
            ind_temp = N + (1:N_ii);
            ind_list(ind_temp) = repmat(ind, M_ii, 1); % spatial index
            a_list(ind_temp) = reshape(repmat(M+(1:M_ii), nxy_data, 1), N_ii, 1); % input index
            v_list(ind_temp) = reshape(data, N_ii, 1);
            N = N + N_ii; % number of non-zero elements in matrix B up to the ii-th part
            M = M + M_ii; % number of columns in matrix B up to the ii-th part
        else
            % If the rectangle of B_struct.pos is a single vertical slice or if it spans the full height of syst.epsilon, we can stack reshaped B_struct(ii).data with zeros to avoid the overhead of building B with index-value pairs. But then B must be built incrementally, which is slow when numel(B_struct) is large.
            nxy_before = sub2ind([ny, nx], m1, n1) - 1;
            nxy_after  = nxy - sub2ind([ny, nx], m2, n2);
            B = [B, [sparse(nxy_before, M_ii); reshape(data, nxy_data, M_ii); sparse(nxy_after, M_ii)]];
        end
    end
    if use_iv_pairs
        B = sparse(ind_list, a_list, v_list, nxy, M);
        if opts.clear_BC; clear data B_struct m_list n_list ind ind_list a_list v_list ind_temp; end
    elseif opts.clear_BC
        clear data B_struct
    end
end
if opts.clear_BC
    B_name = inputname(2); % name of the variable we call B in the caller's workspace; will be empty if there's no variable for it in the caller's workspace
    if ~isempty(B_name)
        evalin('caller', ['clear ', B_name]); % do 'clear B' in caller's workspace
    end
end

% Build the output matrix C from its non-zero elements specified by user
if isstruct(C)
    C_struct = C;
    if ~isfield(C_struct, 'data')
        error('Input argument ''C'' must have field ''data'' when C is a structure array.');
    elseif ~isfield(C_struct, 'pos') && ~isfield(C_struct, 'ind')
        error('Input argument ''C'' must have field ''pos'' or ''ind'' when C is a structure array.');
    elseif isfield(C_struct, 'pos') && isfield(C_struct, 'ind')
        error('Input argument ''C'' cannot have both field ''pos'' and field ''ind'' when C is a structure array.');
    end
    if isfield(C_struct, 'pos')
        for ii = 1:numel(C_struct)
            pos = C_struct(ii).pos;
            if ~(isreal(pos) && isvector(pos) && numel(pos)==4 && isequal(pos, round(pos)) && min(pos)>0)
                error('C(%d).pos must be a positive integer vector with 4 elements.', ii);
            end
        end
    end

    % We first pick the most efficient way to build matrix C.
    % If all of the following are satisfied: (1) numel(C_struct) is small, (2) C_struct.pos is used, and (3) the rectangle specified by C_struct.pos is a single vertical slice or if it spans the full height of syst.epsilon, then we will stack reshaped C_struct(ii).data with zeros. This avoids the overhead of building C with index-value pairs.
    % If any of the above is not satisfied, we will build C with index-value pairs.
    use_iv_pairs = false;
    if numel(C_struct) > 10
        use_iv_pairs = true;
    elseif ~isfield(C_struct, 'pos')
        use_iv_pairs = true;
    else
        for ii = 1:numel(C_struct)
            pos = C_struct(ii).pos;
            if ~(pos(4) == 1 || pos(3) == ny)
                use_iv_pairs = true;
            end
        end
    end

    if use_iv_pairs
        % Construct matrix C from the complete set of index-value pairs
        N_tot = 0; % total number of non-zero elements in C
        for ii = 1:numel(C_struct)
            N_tot = N_tot + numel(C_struct(ii).data);
        end
        ind_list = zeros(N_tot, 0);
        a_list = zeros(N_tot, 0);
        v_list = zeros(N_tot, 0);
        N = 0;
        M = 0;
    else
        % Build matrix C incrementally
        C = sparse(0, nxy);
    end

    % Loop over the structure array
    for ii = 1:numel(C_struct)
        data = C_struct(ii).data;
        if isfield(C_struct, 'pos')
            % C_struct(ii).pos specifies a rectangle inside syst.epsilon; (m1,n1) and (m2,n2) are its two diagonal corners
            pos = C_struct(ii).pos;
            m1 = pos(1); % first index in y
            n1 = pos(2); % first index in x
            m2 = m1 + pos(3) - 1; % last index in y
            n2 = n1 + pos(4) - 1; % last index in x
            nxy_data = pos(3)*pos(4); % number of elements in this rectangle
            if m1 > ny
                error('C(%d).pos(1) = %d exceeds ny = size(syst.epsilon, 1) = %d.', ii, m1, ny);
            elseif n1 > nx
                error('C(%d).pos(2) = %d exceeds nx = size(syst.epsilon, 2) = %d.', ii, n1, nx);
            elseif m2 > ny
                error('C(%d).pos(1) + C(%d).pos(3) - 1 = %d exceeds ny = size(syst.epsilon, 1) = %d.', ii, ii, m2, ny);
            elseif n2 > nx
                error('C(%d).pos(2) + C(%d).pos(4) - 1 = %d exceeds nx = size(syst.epsilon, 2) = %d.', ii, ii, n2, nx);
            elseif ~(ndims(data) <= 3 && isnumeric(data))
                error('C(%d).data must be a 2D or 3D numeric array when C.pos is given.', ii);
            elseif isequal(size(data, [1,2]), [pos(3), pos(4)]) % this way to use size() is supported starting MATLAB R2019b
                M_ii = size(data, 3); % number of outputs
            elseif size(data, 1) == nxy_data
                M_ii = size(data, 2); % number of outputs
            else
                error('size(C(%d).data) = [%d, %d, %d] is not compatible with [C(%d).pos(3), C(%d).pos(4)] = [%d, %d].', ii, size(data, 1), size(data, 2), size(data, 3), ii, ii, pos(3), pos(4));
            end
            if use_iv_pairs
                % convert to linear indices
                m_list = repmat((m1:m2).', 1, pos(4));
                n_list = repmat( n1:n2   , pos(3), 1);
                ind = reshape(sub2ind([ny, nx], m_list, n_list), nxy_data, 1);
            end
        else
            % C_struct(ii).ind specifies linear indices of syst.epsilon
            ind = C_struct(ii).ind;
            if ~(isreal(ind) && isvector(ind) && isequal(ind, round(ind)) && min(ind)>0 && max(ind)<=nxy && numel(ind)==numel(unique(ind)))
                error('C(%d).ind must be a vector with no repeated elements, where every element is an integer between 1 and nx*ny = %d.', ii, nxy);
            elseif ~(ismatrix(data) && isnumeric(data))
                error('C(%d).data must be a 2D numeric matrix when C.ind is given.', ii);
            elseif size(data, 1) ~= numel(ind)
                error('size(C(%d).data, 1) = %d does not match numel(C(%d).ind) = %d', ii, size(data, 1), ii, numel(ind));
            elseif ~iscolumn(ind)
                ind = ind(:);
            end
            nxy_data = numel(ind); % number of linear indices
            M_ii = size(data, 2); % number of outputs
        end

        if use_iv_pairs
            % Build index-value pairs: (a_list, ind_list, v_list)
            N_ii = nxy_data*M_ii; % number of non-zero elements in the ii-th part of matrix C
            ind_temp = N + (1:N_ii);
            ind_list(ind_temp) = repmat(ind, M_ii, 1); % spatial index
            a_list(ind_temp) = reshape(repmat(M+(1:M_ii), nxy_data, 1), N_ii, 1); % input index
            v_list(ind_temp) = reshape(data, N_ii, 1);
            N = N + N_ii; % number of non-zero elements in matrix C up to the ii-th part
            M = M + M_ii; % number of rows in matrix C up to the ii-th part
        else
            % If the rectangle of C_struct.pos is a single vertical slice or if it spans the full height of syst.epsilon, we can stack reshaped C_struct(ii).data with zeros to avoid the overhead of building C with index-value pairs. But then C must be built incrementally, which is slow when numel(C_struct) is large.
            nxy_before = sub2ind([ny, nx], m1, n1) - 1;
            nxy_after  = nxy - sub2ind([ny, nx], m2, n2);
            C = [C; [sparse(M_ii, nxy_before), reshape(data, nxy_data, M_ii).', sparse(M_ii, nxy_after)]];
        end
    end
    if use_iv_pairs
        C = sparse(a_list, ind_list, v_list, M, nxy);
        if opts.clear_BC; clear data C_struct m_list n_list ind ind_list a_list v_list ind_temp; end
    elseif opts.clear_BC
        clear data C_struct
    end
end
if ~isempty(C) && opts.clear_BC
    C_name = inputname(3); % name of the variable we call C in the caller's workspace; will be empty if there's no variable for it in the caller's workspace
    if ~isempty(C_name)
        evalin('caller', ['clear ', C_name]); % do 'clear C' in caller's workspace
    end
end

% Check matrix sizes
[sz_B_1, sz_B_2] = size(B);
if sz_B_1~=nxy; error('size(B,1) must equal nx*ny; size(B,1) = %d, nx*ny = %d.', sz_B_1, nxy); end
if ~opts.return_field_profile && ~use_transpose_B
    [sz_C_1, sz_C_2] = size(C);
    if sz_C_2~=nxy; error('size(C,2) must equal nx*ny; size(C,2) = %d, nx*ny = %d.', sz_C_2, nxy); end
elseif ~isempty(D)
    sz_C_1 = sz_B_2; % to-be used before subtracting D, taking C = transpose(B)
end

t2 = clock; timing_build_BC = etime(t2,t1);
if opts.verbal; fprintf('elapsed time: %7.3f secs\nBuilding A  ... ', timing_build_BC); end

%% Part 2.2: Build matrix A by calling mesti_build_fdfd_matrix()

if (sz_B_2 == 0 || (~opts.return_field_profile && ~use_transpose_B && sz_C_1 == 0)) && ~(isfield(opts, 'store_ordering') && opts.store_ordering)
    % No need to build A if numel(S) = 0 and we don't need to keep the ordering
    A = sparse(nxy, nxy);
    is_symmetric_A = true;
else
    % Build the finite-difference differential operator
    k0dx = (2*pi/syst.wavelength)*(syst.dx);
    [A, is_symmetric_A, xPML, yPML] = mesti_build_fdfd_matrix(syst.epsilon, k0dx, xBC, yBC, xPML, yPML, use_UPML);

    % Use self-energy (instead of PML) to implement exact outgoing boundary condition.
    % self-energy is the retarded Green's function of the surrounding space (with a Dirichlet boundary surrounding the scattering region) evaluated on the surface.
    % Specifically, self-energy = V_SF*inv(A_F)*V_FS where F denotes the surrounding infinite free space, S denotes the scattering region as given by syst.epsilon, A_F is the differential operator of the free space with a Dirichlet boundary surrounding the scattering region S, and V_FS is the coupling matrix between F and S (which is nonzero only along the interface between F and S).
    % On the sides where self-energy is used, we must have Dirichlet boundary condition with no PML. But it's tricky to check which side(s) self-energy is used, so we don't check it.
    % self_energy should have the same symmetry as A; we don't check it.
    if isfield(syst, 'self_energy') && ~isempty(syst.self_energy)
        if ~issparse(syst.self_energy); error('syst.self_energy must be a sparse matrix, if given.'); end
        [sz_SE_1, sz_SE_2] = size(syst.self_energy);
        if sz_SE_1~=nxy || sz_SE_2~=nxy
            error('syst.self_energy must be a square matrix with size nx*ny, if given; size(syst.self_energy) = (%d, %d), nx*ny = %d.', sz_SE_1, sz_SE_2, nxy);
        end
        A = A - syst.self_energy;
    end
end
if opts.clear_syst
    clear syst
    syst_name = inputname(1); % name of the variable we call syst in the caller's workspace; will be empty if there's no variable for it in the caller's workspace
    if ~isempty(syst_name)
        evalin('caller', ['clear ', syst_name]); % do 'clear syst' in caller's workspace
    end
end
opts.is_symmetric_A = is_symmetric_A;

t3 = clock; timing_build_A = etime(t3,t2);
if opts.verbal; fprintf('elapsed time: %7.3f secs\n', timing_build_A); end

%% Part 3: Compute C*inv(A)*B or inv(A)*B by calling mesti_matrix_solver()

% This is where most of the computation is done
% Note that A, B, C in this workspace may be cleared after calling mesti_matrix_solver() if opts.clear_memory = true (which is the default), so we should no longer use A, B, C below
[S, stat] = mesti_matrix_solver(A, B, C, opts);

stat.timing.init = stat.timing.init + timing_init;
stat.timing.build = stat.timing.build + timing_build_BC + timing_build_A; % combine with build time for K
if ~isempty(xPML); stat.xPML = xPML; end
if ~isempty(yPML); stat.yPML = yPML; end

t1 = clock;

% Include the prefactor
S = (opts.prefactor)*S;

if opts.return_field_profile
    % Reshape each of the sz_B_2 field profiles from a vector to a matrix
    S = reshape(S, [ny, nx, sz_B_2]);
end

% subtract D
if ~isempty(D)
    if opts.return_field_profile; error('Input argument ''D'' must be empty for field profile computations where C = [].'); end
    [sz_D_1, sz_D_2] = size(D);
    if sz_D_1~=sz_C_1; error('size(D,1) must equal size(C,1); size(D,1) = %d, size(C,1) = %d.', sz_D_1, sz_C_1); end
    if sz_D_2~=sz_B_2; error('size(D,2) must equal size(B,2); size(D,2) = %d, size(B,2) = %d.', sz_D_2, sz_B_2); end
    S = S - D;
end

t2 = clock; stat.timing.solve = stat.timing.solve + etime(t2,t1); % Add the little bit of post-processing time

stat.timing.total = etime(t2,t0);
if opts.verbal && ~called_from_mesti2s; fprintf('          Total elapsed time: %7.3f secs\n', stat.timing.total); end

end
