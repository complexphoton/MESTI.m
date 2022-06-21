function [A, is_symmetric_A, xPML, yPML] = mesti_build_fdfd_matrix(eps_or_inv_eps, k0dx, xBC, yBC, xPML, yPML, use_UPML)
%MESTI_BUILD_FDFD_MATRIX The finite-difference frequency-domain operator in 2D.
%   A = mesti_build_fdfd_matrix(epsilon, k0dx, xBC, yBC, xPML, yPML, use_UPML)
%   returns A as a sparse matrix representing wave operator
%      [- (d/dx)^2 - (d/dy)^2 - k0^2*epsilon(x,y)]*(dx^2)
%   for the Ez component of transverse-magnetic (TM) fields (Hx, Hy, Ez),
%   discretized on a square grid with grid size dx through center difference.
%   Matrix A has size [ny_Ez*nx_Ez, ny_Ez*nx_Ez].
%
%   A = mesti_build_fdfd_matrix(inv_epsilon, k0dx, xBC, yBC, xPML, yPML, use_UPML)
%   returns A as a sparse matrix representing wave operator
%      [- (d/dx)*(1/epsilon(x,y))_yy*(d/dx) - (d/dy)*(1/epsilon(x,y))_xx*(d/dy)
%       + (d/dy)*(1/epsilon(x,y))_xy*(d/dx) + (d/dx)*(1/epsilon(x,y))_yx*(d/dy)
%       - k0^2]*(dx^2)
%   for the Hz component of transverse-electric (TE) fields (Ex, Ey, Hz). Matrix
%   A has size [ny_Hz*nx_Hz, ny_Hz*nx_Hz].
%
%   === Input Arguments ===
%   eps_or_inv_eps (matrix or cell array; required):
%      When eps_or_inv_eps is a matrix, TM polarization is considered, and
%      epsilon = eps_or_inv_eps is a ny_Ez-by-nx_Ez numeric matrix (real or
%      complex) discretizing the relative permittivity profile epsilon(x,y).
%      Specifically, epsilon(m,n) is the scalar epsilon(x,y) averaged over a
%      square with area dx^2 centered at the point (x_n, y_m) where Ez(x,y) is
%      located on the Yee lattice. It is the zz component of the discretized
%      epsilon(x,y) tensor from subpixel smoothing, used by TM fields. We choose
%         (x_n, y_m) = (n-0.5, m-0.5)*dx,
%         with n = 1, ..., nx_Ez, m = 1, ..., ny_Ez.
%      such that the lower corner of the first pixel epsilon(m=1,n=1) is at
%      (x,y) = (0,0).
%         When eps_or_inv_eps is a cell array, TE polarization is considered,
%      and inv_epsilon = eps_or_inv_eps with
%         inv_epsilon{1}: (1/epsilon(x,y))_xx, size [ny_Ez, nx_Hz]
%         inv_epsilon{2}: (1/epsilon(x,y))_yy, size [ny_Hz, nx_Ez]
%         inv_epsilon{3}: (1/epsilon(x,y))_xy, size [ny_Ez, nx_Ez]
%      The third element, inv_epsilon{3}, is optional and is treated as zero
%      when inv_epsilon only has two elements. The yx component is not
%      specified since we only consider symmetric 1/epsilon tensors where
%      (1/epsilon)_yx = (1/epsilon)_xy.
%         The different components are located at different points:
%         - Hz(x,y) at (x_{n-0.5}, y_{m-0.5}); size [ny_Hz, nx_Hz].
%         - (1/epsilon(x,y))_xx and Dx ~ dHz/dy at (x_{n+0.5}, y_m).
%         - (1/epsilon(x,y))_yy and Dy ~ dHz/dx at (x_n, y_{m+0.5}).
%         - (1/epsilon(x,y))_xy at (x_n, y_m), same as Ez.
%      Here, (x_n, y_m) is the location of Ez and epsilon above.
%         inv_epsilon{1}, inv_epsilon{2}, and inv_epsilon{3} should each be
%      (1/epsilon(x,y))_xx, (1/epsilon(x,y))_yy, and (1/epsilon(x,y))_xy from
%      subpixel smoothing, averaged over a square with area dx^2
%      centered at the points where each of them is located. The
%      1/epsilon(x,y) tensor from subpixel smoothing is given by Eq. (1) of
%      Farjadpour et al, Optics Letters 31, 2972 (2006). Where these points
%      start and end depend on the boundary condition.
%         For periodic, Bloch periodic, and PMCPEC boundary conditions in x,
%      nx_Hz = nx_Ez, and all of the sites on x_{n+0.5} are half a pixel
%      after the corresponding sites on x_n.
%         For PEC boundary condition in x, nx_Hz = nx_Ez + 1, and the sites
%      on x_{n+0.5} start from half a pixel before the first site of x_n and
%      end on half a pixel after the last site of x_n.
%         For PMC boundary condition in x, nx_Hz = nx_Ez - 1, and the sites
%      on x_{n+0.5} start from half a pixel after the first site of x_n and
%      end on half a pixel before the last site of x_n.
%         For PECPMC boundary condition in x, nx_Hz = nx_Ez, and all of the
%      sites on x_{n+0.5} are half a pixel before the corresponding sites on
%      x_n.
%         Similar applies to boundary conditions in y.
%   k0dx (numeric scalar, real or complex; required):
%      Normalized frequency k0*dx = (2*pi/vacuum_wavelength)*dx.
%   xBC (character vector or numeric scalar; required):
%      Boundary condition (BC) at the two ends in x direction, effectively
%      specifying Ez(m,n) or Hz(m,n) at n=0 and n=nx_Ez+1 or nx_Hz+1, one
%      pixel beyond the computation domain. For character inputs, the options
%      are:
%         'periodic' - Ez(m,n+nx_Ez) = Ez(m,n)
%                      Hz(m,n+nx_Hz) = Hz(m,n)
%         'PEC'      - Ez(m,0) = Ez(m,nx_Ez+1) = 0
%                      Hz(m,0) = Hz(m,1); Hz(m,nx_Hz+1) = Hz(m,nx_Hz)
%         'PMC'      - Ez(m,0) = Ez(m,1); Ez(m,nx_Ez+1) = Ez(m,nx_Ez)
%                      Hz(m,0) = Hz(m,nx_Hz+1) = 0
%         'PECPMC'   - Ez(m,0) = 0; Ez(m,nx_Ez+1) = Ez(m,nx_Ez)
%                      Hz(m,0) = Hz(m,1); Hz(m,nx_Hz+1) = 0
%         'PMCPEC'   - Ez(m,0) = Ez(m,1); Ez(m,nx_Ez+1) = 0
%                      Hz(m,0) = 0; Hz(m,nx_Hz+1) = Hz(m,nx_Hz)
%      where PEC stands for perfect electric conductor (for which Ez = 0 and
%      Ey ~ dHz/dx = 0 at the boundary) and PMC stands for perfect magnetic
%      conductor (for which Hz = 0 and Hy ~ dEz/dx = 0 at the boundary).
%         When xBC is a numeric scalar, the Bloch periodic boundary condition is
%      used with
%                    - Ez(m,n+nx_Ez) = Ez(m,n)*exp(1i*x_BC)
%                      Hz(m,n+nx_Hz) = Hz(m,n)*exp(1i*x_BC)
%      In other words, xBC = kx_B*p where kx_B is the Bloch wave number, and p =
%      nx_Ez*dx or nx_Hz*dx is the periodicity in x.
%   yBC (character vector or numeric scalar; required):
%      Boundary condition in y direction, analogous to xBC.
%   xPML (two-element cell array or scalar structure or []; optional):
%      Parameters for perfectly matched layer (PML) in x direction.
%      xPML = [] or {[], []} means no PML on either side (default).
%      xPML = {[], PML} means no PML on the left.
%      xPML = {PML, []} means no PML on the right.
%      xPML = {PML1, PML2} means PML on both sides.
%      xPML = PML means PML on both sides with the same parameters, equivalent
%         to xPML = {PML, PML}.
%      In each case, PML is a scalar structure with the following fields:
%         npixels (non-negative integer scalar; required): Number of PML pixels.
%            Note this is within epsilon or inv_epsilon.
%         power_sigma (non-negative scalar; optional): Power of the polynomial 
%            grading for the conductivity sigma; defaults to 3.
%         sigma_max_over_omega (non-negative scalar; optional):
%            Conductivity at the end of the PML; defaults to
%               0.8*(power_sigma+1)/(k0dx*sqrt(epsilon_bg)).
%            where epsilon_bg is the average relative permittivity along the
%            last slice of the PML. This is used to attenuate propagating waves.
%         power_kappa (non-negative scalar; optional): Power of the polynomial
%            grading for the real-coordinate-stretching factor kappa; defaults
%            to 3.
%         kappa_max (real scalar no smaller than 1; optional):
%            Real-coordinate-stretching factor at the end of the PML; defaults
%            to 15. This is used to accelerate the attenuation of evanescent 
%            waves. kappa_max = 1 means no real-coordinate stretching.
%         power_alpha (non-negative scalar; optional): Power of the polynomial
%            grading for the CFS alpha factor; defaults to 1.
%         alpha_max_over_omega (non-negative scalar; optional): Complex-
%            frequency-shifting (CFS) factor at the beginning of the PML. This
%            is typically used in time-domain simulations to suppress late-time
%            (low-frequency) reflections. We don't use it by default 
%            (alpha_max_over_omega = 0) since we are in frequency domain.
%      We use the following PML coordinate-stretching factor:
%         s(u) = kappa(u) + sigma(u)./(alpha(u) - i*omega)
%      with
%         sigma(u)/omega = sigma_max_over_omega*(u.^power_sigma),
%         kappa(u) = 1 + (kappa_max-1)*(u.^power_kappa),
%         alpha(u)/omega = alpha_max_over_omega*((1-u).^power_alpha),
%      where omega is frequency, and u goes linearly from 0 at the beginning of
%      the PML to 1 at the end of the PML. 
%   yPML (two-element cell array or scalar structure or []; optional):
%      Parameters for PML in y direction, analogous to xPML.
%   use_UPML (logical scalar; optional, defaults to true):
%      Whether to use uniaxial PML (UPML) or not. If not, stretched-coordinate
%      PML (SC-PML) will be used.
%
%   === Output Arguments ===
%   A (sparse matrix):
%      Sparse matrix discretizing the differential operator on Ez or Hz. It has
%      size [ny_Ez*nx_Ez, ny_Ez*nx_Ez] for TM, size [ny_Hz*nx_Hz, ny_Hz*nx_Hz]
%      for TE.
%   is_symmetric_A (logical scalar):
%      Whether matrix A is symmetric or not.
%   xPML (two-element cell array or []):
%      PML parameters used on the low and high sides of x direction, if any.
%   yPML (two-element cell array or []):
%      PML parameters used on the low and high sides of y direction, if any.

% Check input arguments
% xBC and yBC are checked in build_ddx_E(); xPML and yPML are checked in set_PML_params()
if nargin < 4; error('Not enough input arguments.'); end
if ~iscell(eps_or_inv_eps)
    use_TM = true;
    epsilon = eps_or_inv_eps;
    if ~(ismatrix(epsilon) && isnumeric(epsilon))
        error('Input argument epsilon must be a numeric matrix.');
    end
else
    use_TM = false;
    inv_epsilon = eps_or_inv_eps;
    if numel(inv_epsilon) ~= 2 && numel(inv_epsilon) ~= 3
        error('Input argument inv_epsilon must be a two-element or three-element cell array.');
    elseif ~(ismatrix(inv_epsilon{1}) && isnumeric(inv_epsilon{1}))
        error('inv_epsilon{1} must be a numeric matrix.');
    elseif ~(ismatrix(inv_epsilon{2}) && isnumeric(inv_epsilon{2}))
        error('inv_epsilon{2} must be a numeric matrix.');
    elseif numel(inv_epsilon) == 3 && ~(ismatrix(inv_epsilon{3}) && isnumeric(inv_epsilon{3}))
        error('inv_epsilon{3} must be a numeric matrix, if given.');
    end
end
if ~(isscalar(k0dx) && isnumeric(k0dx)); error('Input argument k0dx must be a numeric scalar.'); end

% Assign default value for optional inputs
if nargin < 5; xPML = []; end
if nargin < 6; yPML = []; end
if nargin < 7; use_UPML = true; end % use UPML instead of SC-PML by default
if ~(islogical(use_UPML) && isscalar(use_UPML))
    error('Input argument use_UPML must be a logical scalar, if given.');
end

% Number of sites in y and x
if use_TM
    [ny_Ez, nx_Ez] = size(epsilon);
    nxy = nx_Ez*ny_Ez;
else
    [ny_Ez, nx_Hz_temp] = size(inv_epsilon{1}); % (1/epsilon)_xx
    [ny_Hz_temp, nx_Ez] = size(inv_epsilon{2}); % (1/epsilon)_yy
    nxy = nx_Hz_temp*ny_Hz_temp;
end

% Nothing to do if the matrix size is zero
if nxy == 0
    A = sparse(nxy, nxy);
    is_symmetric_A = true;
    return
end

% Estimate background permittivity, used to assign default sigma_max_over_omega for PML
epsilon_bg_x = [1, 1];
epsilon_bg_y = [1, 1];
if ~isempty(xPML)
    if use_TM
        epsilon_bg_x = real([mean(epsilon(:,1)), mean(epsilon(:,end))]);

        % Make sure that the two sides are the same when the system is periodic; this ensures that the s-factor is continuous across the periodic boundary.
        if isnumeric(xBC) || strcmpi(xBC, 'periodic')
            epsilon_bg_x = mean(epsilon_bg_x)*[1, 1];
        end
    else
        % If nx_Ez is zero (possible when nx_Hz = 1 with PEC boundary), we can't take the mean, but in such case there must be no PML in x.
        if nx_Ez ~= 0
            epsilon_bg_x = real(1./[mean(inv_epsilon{2}(:,1)), mean(inv_epsilon{2}(:,end))]);
        end

        % Make sure that the two sides are the same when the system is periodic; this ensures that the s-factor is continuous across the periodic boundary.
        % Under periodic or Bloch BC for TE fields, Ey and (1/epsilon)_yy are placed half a pixel in x before Hz, so the first slice of (1/epsilon)_yy is half on the left, half on the right, and already contains spatial average of the two sides.
        if isnumeric(xBC) || strcmpi(xBC, 'periodic')
            epsilon_bg_x = epsilon_bg_x(1)*[1, 1];
        end
    end
end
if ~isempty(yPML)
    if use_TM
        epsilon_bg_y = real([mean(epsilon(1,:)), mean(epsilon(end,:))]);
        if isnumeric(yBC) || strcmpi(yBC, 'periodic')
            epsilon_bg_y = mean(epsilon_bg_y)*[1, 1];
        end
    else
        if ny_Ez ~= 0
            epsilon_bg_y = real(1./[mean(inv_epsilon{1}(1,:)), mean(inv_epsilon{1}(end,:))]);
        end
        if isnumeric(yBC) || strcmpi(yBC, 'periodic')
            epsilon_bg_y = epsilon_bg_y(1)*[1, 1];
        end
    end
end

% Set default values for PML parameters
xPML = set_PML_params(xPML, k0dx, epsilon_bg_x, 'x');
yPML = set_PML_params(yPML, k0dx, epsilon_bg_y, 'y');

% Build the first derivative and average matrices on E, such that
%   ddx_E*f = (df/dx)*dx
%   avg_x_E*f = average of f among two neighboring pixels
% where f = Ey or Ez is a 1d vector.
% Note that sx_E and sx_H are column vectors.
[ddx_E, avg_x_E, sx_E, sx_H] = build_ddx_E(nx_Ez, xBC, xPML, 'x', use_TM); % ddx_E operates on Ey or Ez
[ddy_E, avg_y_E, sy_E, sy_H] = build_ddx_E(ny_Ez, yBC, yPML, 'y', use_TM); % ddy_E operates on Ex or Ez

% The derivative matrices on H
ddx_H = -ddx_E'; % ddx_H operates on Hy or Hz
ddy_H = -ddy_E'; % ddy_H operates on Hx or Hz

% Number of sites of Ez
% The two lines below are redundant; they are kept for clarity.
nx_Ez = size(ddx_E, 2); % nx_Ez = nx_Ey = nx_Hx
ny_Ez = size(ddy_E, 2); % ny_Ez = ny_Ex = ny_Hy

% Number of sites of Hz
nx_Hz = size(ddx_E, 1); % nx_Hz = nx_Hy = nx_Ex
ny_Hz = size(ddy_E, 1); % ny_Hz = ny_Hx = ny_Ey

if use_TM
    % Include coordinate stretching into the first derivatives
    ddx_E = spdiag(1./sx_H, nx_Hz)*ddx_E;
    ddy_E = spdiag(1./sy_H, ny_Hz)*ddy_E;
    if ~use_UPML
        ddx_H = spdiag(1./sx_E, nx_Ez)*ddx_H;
        ddy_H = spdiag(1./sy_E, ny_Ez)*ddy_H;
    end
else
    % Include coordinate stretching into the first derivatives
    ddx_H = spdiag(1./sx_E, nx_Ez)*ddx_H;
    ddy_H = spdiag(1./sy_E, ny_Ez)*ddy_H;
    if ~use_UPML
        ddx_E = spdiag(1./sx_H, nx_Hz)*ddx_E;
        ddy_E = spdiag(1./sy_H, ny_Hz)*ddy_E;
    end

    % Check if size(inv_epsilon{1}) and size(inv_epsilon{2}) are compatible with the boundary conditions
    if nx_Hz_temp ~= nx_Hz
        error('nx_Hz = size(inv_epsilon{1},2) = %d but should be %d instead given nx_Ez = size(inv_epsilon{2},2) = %d and the boundary condition in x.', nx_Hz_temp, nx_Hz, nx_Ez);
    elseif ny_Hz_temp ~= ny_Hz
        error('ny_Hz = size(inv_epsilon{2},1) = %d but should be %d instead given ny_Ez = size(inv_epsilon{1},1) = %d and the boundary condition in y.', ny_Hz_temp, ny_Hz, ny_Ez);
    end
    inv_epsilon_xx = spdiag(inv_epsilon{1}(:), nx_Hz*ny_Ez);
    inv_epsilon_yy = spdiag(inv_epsilon{2}(:), nx_Ez*ny_Hz);

    % off-diagonal components of 1/epsilon
    if numel(inv_epsilon) == 3 && nnz(inv_epsilon{3}) > 0
        include_inv_epsilon_xy = true;

        % Check if size(inv_epsilon{3}) is compatible with size(inv_epsilon{1}) and size(inv_epsilon{2})
        [ny_Ez_temp, nx_Ez_temp] = size(inv_epsilon{3}); % inv_epsilon_xy
        if nx_Ez_temp ~= nx_Ez
            error('nx_Ez = size(inv_epsilon{3},2) = %d does not match nx_Ez = size(inv_epsilon{2},2) = %d.', nx_Ez_temp, nx_Ez);
        elseif ny_Ez_temp ~= ny_Ez
            error('ny_Ez = size(inv_epsilon{3},1) = %d does not match ny_Ez = size(inv_epsilon{1},1) = %d.', ny_Ez_temp, ny_Ez);
        end
        inv_epsilon_xy = spdiag(inv_epsilon{3}(:), nx_Ez*ny_Ez);
        inv_epsilon_yx = inv_epsilon_xy; % 1/epsilon is a symmetric tensor

        % Matrix performing average on Hz (or Dx, Dy)
        avg_x_H = avg_x_E';
        avg_y_H = avg_y_E';
    else
        include_inv_epsilon_xy = false;
    end
end

% Build the operators; use Kronecker outer product to go from 1D to 2D.
% TM: A = [- (d/dx)^2 - (d/dy)^2 - k0^2*epsilon(x,y)]*(dx^2)
% TE: A = [- (d/dx)*(1/epsilon(x,y))_yy*(d/dx) - (d/dy)*(1/epsilon(x,y))_xx*(d/dy)
%          + (d/dy)*(1/epsilon(x,y))_xy*(d/dx) + (d/dx)*(1/epsilon(x,y))_yx*(d/dy) - k0^2]*(dx^2)
if ~use_UPML
    % Stretched-coordinate PML (SC-PML); here, both 1/s_E and 1/s_H have already been multiplied onto the first derivatives.
    if use_TM
        A = - kron(ddx_H*ddx_E, speye(ny_Ez)) ...
            - kron(speye(nx_Ez), ddy_H*ddy_E) ...
            - spdiag((k0dx^2)*epsilon(:), nxy);
    else
        A = - kron(ddx_E, speye(ny_Hz)) * inv_epsilon_yy * kron(ddx_H, speye(ny_Hz)) ...
            - kron(speye(nx_Hz), ddy_E) * inv_epsilon_xx * kron(speye(nx_Hz), ddy_H) ...
            - (k0dx^2)*speye(nxy);
        if include_inv_epsilon_xy
            % Following Oskooi et al, Optics Letters 34, 2778 (2009), we average two points of Dy ~ dHz/dx along y, multiply by (1/epsilon)_xy to get Ex, and then average two such points along x to get Ex on the site of Yee lattice we need before taking Hz ~ dEx/dy. The average in y and derivative in x commute. Similarly for the other term. To summarize:
            %   first term: (d/dy)*avg_x*(1/epsilon)_xy*avg_y*(d/dx)*Hz
            %  second term: (d/dx)*avg_y*(1/epsilon)_yx*avg_x*(d/dy)*Hz
            % Many of these matrix elements cancel among the two terms, so we sum them first before adding to A.
            A = A + (kron(avg_x_E, ddy_E) * inv_epsilon_xy * kron(ddx_H, avg_y_H) ...
                   + kron(ddx_E, avg_y_E) * inv_epsilon_yx * kron(avg_x_H, ddy_H));
        end
    end
else
    % Uniaxial PML (UPML); here, the first 1/s has already been multiplied onto the first derivatives.
    % A_UPML = (sx_E*sy_E)*A_SCPML for TM fields
    % A_UPML = (sx_H*sy_H)*A_SCPML for TE fields
    if use_TM
        sxy_E = sy_E.*reshape(sx_E, 1, nx_Ez); % use implicit expansion
        A = - kron(ddx_H*ddx_E, spdiag(sy_E, ny_Ez)) ...
            - kron(spdiag(sx_E, nx_Ez), ddy_H*ddy_E) ...
            - spdiag((k0dx^2)*(sxy_E(:).*epsilon(:)), nxy);
    else
        sxy_H = sy_H.*reshape(sx_H, 1, nx_Hz); % use implicit expansion
        sx_H = spdiag(sx_H, nx_Hz);
        sy_H = spdiag(sy_H, ny_Hz);
        A = - kron(ddx_E, sy_H) * inv_epsilon_yy * kron(ddx_H, speye(ny_Hz)) ...
            - kron(sx_H, ddy_E) * inv_epsilon_xx * kron(speye(nx_Hz), ddy_H) ...
            - spdiag((k0dx^2)*sxy_H(:), nxy);
        if include_inv_epsilon_xy
            A = A + (kron(sx_H*avg_x_E, ddy_E) * inv_epsilon_xy * kron(ddx_H, avg_y_H) ...
                   + kron(ddx_E, sy_H*avg_y_E) * inv_epsilon_yx * kron(avg_x_H, ddy_H));
        end
    end
end

% determine the symmetry of matrix A, assuming no spatial symmetry in epsilon or inv_epsilon
if (isnumeric(xBC) && xBC ~= 0 && xBC ~= pi && nx_Ez > 1) || (isnumeric(yBC) && yBC ~= 0 && yBC ~= pi && ny_Ez > 1)
    % Bloch periodic boundary condition with ka != 0 or pi breaks the symmetry of A because its ddx is complex-valued
    is_symmetric_A = false;
elseif ~isempty(xPML) || ~isempty(yPML)
    if ~use_UPML
        % SC-PML breaks the symmetry of A
        is_symmetric_A = false;
    elseif ~use_TM && include_inv_epsilon_xy
        % off-diagonal components of inv(epsilon) breaks the symmetry of A when PML is used
        is_symmetric_A = false;
    else
        is_symmetric_A = true;
    end
elseif ~use_TM && include_inv_epsilon_xy && ((isnumeric(xBC) && xBC ~= 0 && xBC ~= pi && nx_Hz == 1 && ny_Hz > 1) || (isnumeric(yBC) && yBC ~= 0 && yBC ~= pi && ny_Hz == 1 && nx_Hz > 1))
    % A single layer with Bloch periodic boundary condition can still break the symmetry of A when inv(epsilon) has off-diagonal components
    is_symmetric_A = false;
else
    is_symmetric_A = true;
end

end


% helper function to generate a sparse diagonal matrix
% d must be a column vector
function A = spdiag(d, n)
    A = spdiags(d, 0, n, n);
end


function [ddx, avg, s_E, s_H] = build_ddx_E(n_E, BC, PML, direction, use_TM)
% ddx: first-derivative matrix (acting on Ey or Ez)
% avg: average matrix (acting on Ey or Ez)
% s_E: x-coordinate-stretching factor for Ey or Ez (on integer sites)
% s_H: x-coordinate-stretching factor for Hy or Hz (on half-integer sites)
% n_E: number of sites in x for Ey or Ez

if ~((ischar(BC) && isrow(BC)) || ((isstring(BC) || isnumeric(BC)) && isscalar(BC)))
    error('Input argument %sBC must be a character vector or string, or numeric scalar (for Bloch periodic BC).', direction);
end

% Handle periodic and Bloch periodic boundary conditions
if strcmpi(BC, 'Bloch')
    error('To use Bloch periodic boundary condition, set BC to k_B*p where k_B is the Bloch wave number and p is the periodicity.');
elseif isnumeric(BC)
    ka = BC;
    BC = 'Bloch';
    if ~isreal(ka)
        warning('k%s_B*a = %g + 1i*%g is a complex number.', direction, real(ka), imag(ka));
    end
elseif strcmpi(BC, 'periodic')
    ka = 0;
    BC = 'Bloch';
end

% Build the first-derivative matrix and the average matrix on E
% f = [f(1), ..., f(n_E)].' with f = Ey or Ez, on integer sites
% df = (df/dx)*dx, proportional to Hy or Hz, on half-integer sites
% avg_f = average of f between two neighboring sites, on half-integer sites
if strcmpi(BC, 'Bloch')
    % f(n_E+1) = f(1)*exp(1i*ka); f(0) = f(n_E)*exp(-1i*ka)
    % ddx*f = [df(1.5), ..., df(n_E+0.5)].'
    ddx = spdiags([ones(n_E,1),-ones(n_E,1),exp(1i*ka)*ones(n_E,1)], [1,0,1-n_E], n_E, n_E);
    % avg*f = [avg_f(1.5), ..., avg_f(n_E+0.5)].'
    avg = spdiags([ones(n_E,1),ones(n_E,1),exp(1i*ka)*ones(n_E,1)]/2, [1,0,1-n_E], n_E, n_E);
elseif strcmpi(BC, 'PEC') % PEC on both sides
    % f(0) = f(n_E+1) = 0
    % ddx*f = [df(0.5), ..., df(n_E+0.5)].'
    ddx = spdiags([ones(n_E,1),-ones(n_E,1)], [0,-1], n_E+1, n_E);
    % avg*f = [avg_f(0.5), ..., avg_f(n_E+0.5)].'
    avg = spdiags(ones(n_E,2)/2, [0,-1], n_E+1, n_E);
elseif strcmpi(BC, 'PMC') % PMC on both sides
    % f(0) = f(1); f(n_E+1) = f(n_E)
    % ddx*f = [df(1.5), ..., df(n_E-0.5)].'; exclude df(0.5) and df(n_E+0.5) because they are zero
    ddx = spdiags([ones(n_E-1,1),-ones(n_E-1,1)], [1,0], n_E-1, n_E);
    % avg*f = [avg_f(1.5), ..., avg_f(n_E-0.5)].'; exclude avg_f(0.5) and avg_f(n_E+0.5)
    avg = spdiags(ones(n_E-1,2)/2, [1,0], n_E-1, n_E);
elseif strcmpi(BC, 'PECPMC') % PEC on the low side, PMC on the high side
    % f(0) = 0; f(n_E+1) = f(n_E)
    % ddx*f = [df(0.5), ..., df(n_E-0.5)].'; exclude df(n_E+0.5) because it is zero
    ddx = spdiags([ones(n_E,1),-ones(n_E,1)], [0,-1], n_E, n_E);
    % avg*f = [avg_f(0.5), ..., avg_f(n_E-0.5)].'; we exclude avg_f(n_E+0.5)
    avg = spdiags(ones(n_E,2)/2, [0,-1], n_E, n_E);
elseif strcmpi(BC, 'PMCPEC') % PMC on the low side, PEC on the high side
    % f(0) = f(1); f(n_E+1) = 0
    % ddx*f = [df(1.5), ..., df(n_E+0.5)].'; exclude df(0.5) because it is zero
    ddx = spdiags([ones(n_E,1),-ones(n_E,1)], [1,0], n_E, n_E);
    % avg*f = [avg_f(1.5), ..., avg_f(n_E+0.5)].'; exclude avg_f(0.5)
    avg = spdiags(ones(n_E,2)/2, [1,0], n_E, n_E);
else
    error('Input argument %sBC = %s is not a supported option.', direction, BC);
end

n_H = size(ddx,1); % number of sites in x for df/dx (ie, Hy or Hz)

% coordinate-stretching factor in x: s(x) = kappa(x) + sigma(x)/(alpha(u) - i*omega)
s_E = ones(n_E,1); % s-factor for Ey or Ez (on integer sites)
s_H = ones(n_H,1); % s-factor for Hy or Hz (on half-integer sites)

% no coordinate stretching if no PML is specified
if isempty(PML)
    return
end

% number of PML pixels on the low and high sides
npixels = [PML{1}.npixels, PML{2}.npixels];

% Here, we let f = Ez for TM fields; f = Hz for TE fields
% u_PML_1 and ind_PML_1 acts on f and d^2f/dx^2
% u_PML_2 and ind_PML_2 acts on df/dx
if use_TM
    n = n_E;
else
    n = n_H;
end

% Cannot have more PML pixels than the number of pixels
% In our implementation, the conductivity goes to zero at one pixel before PML. So there needs to be at least one more pixel in addition to PML.
if sum(npixels) >= n
    error('Total number of pixels = %d in %s direction must be greater than the number of PML pixels = %d + %d = %d but is not.', n, direction, npixels(1), npixels(2), sum(npixels));
end

% Below, u(x) is a function that goes linearly from u(x)=0 at one site before PML to u(x)=1 at the "end of PML".
% Note that where the "end of PML" is depends on the boundary condition.
% Let index i=0 be one site before PML, so u(i=0)=0.
% Then, i=1 is the first site of PML and i=npixels is the last site of PML we explicitly simulate. But we do not set u(i=npixels)=1.
% For Dirichlet BC, the BC is such that f=0 at i=npixels+1, so we let the end of PML be i=npixels+1, with an effective PML thickness of npixels+1 pixels.
% For Neumann BC, the BC is such that df=0 at i=npixels+0.5, so we let the end of PML be i=npixels+0.5, with an effective PML thickness of npixels+0.5 pixels.
% For periodic and Bloch periodic BC, we let u(x) be symmetric on the two sides of f with u(0.5)=u(n+0.5)=1; note that f(0.5) and f(n+0.5) are the same site (with a possible Bloch phase difference). If the PML parameters on the two sides are the same, the s-factor will be continuous across the periodic boundary.

% Construct u(x) and their corresponding indices for df/dx
if strcmpi(BC, 'Bloch')
    if npixels(1)*npixels(2) == 0
        warning('Bloch periodic boundary condition is used with a single-sided PML in %s direction; transmission through PML will only undergo single-pass attenuation.', direction);
    end
    npixels_effective = [npixels(1)+0.5, npixels(2)+0.5];
    if use_TM
        % f = Ez; ddx*f = [df(1.5), ..., df(n_E+0.5)].'
        % no s-factor for df(0.5) since we only consider df(n_E+0.5)
        u_PML_2 = {((1:npixels(1))-0.5)/npixels_effective(1), ((1:(npixels(2)+1))-0.5)/npixels_effective(2)};
        ind_PML_2 = {fliplr(1:npixels(1)), (n+1)-fliplr(1:(npixels(2)+1))};
    else
        % f = Hz; (-ddx')*f = [df(0.5), ..., df(n_H-0.5)].'
        % no s-factor for df(n_H+0.5) since we only consider df(0.5)
        u_PML_2 = {((1:(npixels(1)+1))-0.5)/npixels_effective(1), ((1:npixels(2))-0.5)/npixels_effective(2)};
        ind_PML_2 = {fliplr(1:(npixels(1)+1)), (n+1)-fliplr(1:npixels(2))};
    end
elseif (use_TM && strcmpi(BC, 'PEC')) || (~use_TM && strcmpi(BC, 'PMC')) % Dirichlet on both sides
    npixels_effective = [npixels(1)+1, npixels(2)+1];
    u_PML_2 = {((1:(npixels(1)+1))-0.5)/npixels_effective(1), ((1:(npixels(2)+1))-0.5)/npixels_effective(2)};
    ind_PML_2 = {fliplr(1:(npixels(1)+1)), (n+2)-fliplr(1:(npixels(2)+1))};
elseif (use_TM && strcmpi(BC, 'PMC')) || (~use_TM && strcmpi(BC, 'PEC')) % Neumann on both sides
    npixels_effective = [npixels(1)+0.5, npixels(2)+0.5];
    u_PML_2 = {((1:npixels(1))-0.5)/npixels_effective(1), ((1:npixels(2))-0.5)/npixels_effective(2)};
    ind_PML_2 = {fliplr(1:npixels(1)), n-fliplr(1:npixels(2))};
elseif (use_TM && strcmpi(BC, 'PECPMC')) || (~use_TM && strcmpi(BC, 'PMCPEC')) % Dirichlet on the low side, Neumann on the high side
    npixels_effective = [npixels(1)+1, npixels(2)+0.5];
    u_PML_2 = {((1:(npixels(1)+1))-0.5)/npixels_effective(1), ((1:npixels(2))-0.5)/npixels_effective(2)};
    ind_PML_2 = {fliplr(1:(npixels(1)+1)), (n+1)-fliplr(1:npixels(2))};
elseif (use_TM && strcmpi(BC, 'PMCPEC')) || (~use_TM && strcmpi(BC, 'PECPMC')) % Neumann on the low side, Dirichlet on the high side
    npixels_effective = [npixels(1)+0.5, npixels(2)+1];
    u_PML_2 = {((1:npixels(1))-0.5)/npixels_effective(1), ((1:(npixels(2)+1))-0.5)/npixels_effective(2)};
    ind_PML_2 = {fliplr(1:npixels(1)), (n+1)-fliplr(1:(npixels(2)+1))};
else
    error('Input argument %sBC = ''%s'' is not a supported option.', direction, BC);
end

% Construct u(x) and their corresponding indices for f and d^2f/dx^2
u_PML_1 = {(1:npixels(1))/npixels_effective(1), (1:npixels(2))/npixels_effective(2)};
ind_PML_1 = {fliplr(1:npixels(1)), (n+1)-fliplr(1:npixels(2))};

% Recall that f = Ez for TM fields; f = Hz for TE fields
if use_TM
    u_PML_E = u_PML_1;
    ind_PML_E = ind_PML_1;
    u_PML_H = u_PML_2;
    ind_PML_H = ind_PML_2;
else
    u_PML_H = u_PML_1;
    ind_PML_H = ind_PML_1;
    u_PML_E = u_PML_2;
    ind_PML_E = ind_PML_2;
end

% Loop over PML on the two sides
for ii = 1:2
    if npixels(ii) > 0
        % coordinate-stretching factor s(u) from Eq 7.73 of Taflove & Hagness's 2005 FDTD book
        % sigma is the conductivity, equivalent to imag-coordinate stretching, used to attenuate propagating waves.
        % kappa is real-coordinate stretching, used to accelerate the attenuation of evanescent waves.
        % alpha is used for complex frequency shifting (CFS) to suppress reflection of low-frequency components for time-domain simulations.
        % In general, kappa, sigma, alpha can all be arbitrary functions of position.
        % To minimize discretization-induced reflection, kappa should start from 1, and sigma should start from 0.
        % We use a polynomial grading for all of them: Eqs 7.60 and 7.79 of Taflove & Hagness's 2005 FDTD book.
        kappa = @(u) 1 + (PML{ii}.kappa_max-1)*(u.^(PML{ii}.power_kappa));
        sigma_over_omega = @(u) PML{ii}.sigma_max_over_omega*(u.^(PML{ii}.power_sigma));
        alpha_over_omega = @(u) PML{ii}.alpha_max_over_omega*((1-u).^(PML{ii}.power_alpha));
        func_s = @(u) kappa(u) + sigma_over_omega(u)./(alpha_over_omega(u) - 1i);

        % Evaluate s(u) on integer and half-integer sites
        % Note s_E and s_H should be sampled directly from the polynomial functions; if we use the polynomial for one and linear interpolation for the other, the performance will be much worse.
        s_E(ind_PML_E{ii}) = func_s(u_PML_E{ii}).'; % column vector
        s_H(ind_PML_H{ii}) = func_s(u_PML_H{ii}).'; % column vector
    end
end

end


function PML = set_PML_params(PML, k0dx, epsilon_bg, direction)
% Set default values for PML parameters

if isempty(PML)
    return
elseif ~((isstruct(PML) && isscalar(PML)) || (iscell(PML) && numel(PML) == 2))
    error('Input argument %sPML must be a scalar structure or a two-element cell array.', direction)
elseif isstruct(PML)
    % Use the same PML parameters on both sides if only one set of parameters is given
    PML = {PML, PML};
end

% Loop over PML parameters on the two sides
for ii = 1:2
    % Empty means no PML, and we let npixels = 0
    if isempty(PML{ii})
        PML{ii}.npixels = 0;
        continue
    elseif ~(isstruct(PML{ii}) && isscalar(PML{ii}))
        error('%sPML{%d} must be a scalar structure.', direction, ii);
    end

    % Number of PML pixels
    if ~isfield(PML{ii}, 'npixels') || isempty(PML{ii}.npixels)
        error('%sPML{%d} must contain field ''npixels''.', direction, ii);
    else
        temp = PML{ii}.npixels;
        if ~(isreal(temp) && isscalar(temp) && temp >= 0 && round(temp) == temp)
            error('%sPML{%d}.npixels must be a non-negative integer scalar.', direction, ii);
        end
        if temp == 0
            continue
        end
    end

    % Power of polynomial grading for the conductivity sigma
    if ~isfield(PML{ii}, 'power_sigma') || isempty(PML{ii}.power_sigma)
        % From Table 7.1 of Taflove & Hagness's 2005 FDTD book
        PML{ii}.power_sigma = 3;
    else
        temp = PML{ii}.power_sigma;
        if ~(isreal(temp) && isscalar(temp) && temp >= 0)
            error('%sPML{%d}.power_sigma must be a non-negative scalar.', direction, ii);
        end
    end

    % Conductivity at the end of the PML
    if ~isfield(PML{ii}, 'sigma_max_over_omega') || isempty(PML{ii}.sigma_max_over_omega)
        % Eq 7.67 of Taflove & Hagness's 2005 FDTD book
        PML{ii}.sigma_max_over_omega = 0.8*(PML{ii}.power_sigma+1)/(k0dx*sqrt(epsilon_bg(ii)));
    else
        temp = PML{ii}.sigma_max_over_omega;
        if ~(isreal(temp) && isscalar(temp) && temp >= 0)
            error('%sPML{%d}.sigma_max_over_omega must be a non-negative scalar.', direction, ii);
        end
    end

    % Maximal coordinate stretching factor
    if ~isfield(PML{ii}, 'kappa_max') || isempty(PML{ii}.kappa_max)
        % From Table 7.1 of Taflove & Hagness's 2005 FDTD book
        PML{ii}.kappa_max = 15;
    else
        temp = PML{ii}.kappa_max;
        if ~(isreal(temp) && isscalar(temp) && temp >= 1)
            error('%sPML{%d}.kappa_max must be a real scalar that equals or is larger than 1.', direction, ii);
        end
    end

    % Power of polynomial grading for the coordinate stretching factor kappa
    if ~isfield(PML{ii}, 'power_kappa') || isempty(PML{ii}.power_kappa)
        % From Table 7.1 of Taflove & Hagness's 2005 FDTD book
        PML{ii}.power_kappa = 3;
    else
        temp = PML{ii}.power_kappa;
        if ~(isreal(temp) && isscalar(temp) && temp >= 0)
            error('%sPML{%d}.power_kappa must be a non-negative scalar.', direction, ii);
        end
    end

    % Maximal alpha factor for complex frequency shifting (CFS)
    if ~isfield(PML{ii}, 'alpha_max_over_omega') || isempty(PML{ii}.alpha_max_over_omega)
        % CFS is meant to suppress reflection of low-frequency components for time-domain simulations; it is not necessary for frequency-domain simulations
        PML{ii}.alpha_max_over_omega = 0;
    else
        temp = PML{ii}.alpha_max_over_omega;
        if ~(isreal(temp) && isscalar(temp) && temp >= 0)
            error('%sPML{%d}.alpha_max_over_omega must be a non-negative scalar.', direction, ii);
        end
    end

    % Power of polynomial grading for the alpha factor
    if ~isfield(PML{ii}, 'power_alpha') || isempty(PML{ii}.power_alpha)
        % not relevant when alpha_max_over_omega = 0
        PML{ii}.power_alpha = 1;
    else
        temp = PML{ii}.power_alpha;
        if ~(isreal(temp) && isscalar(temp) && temp >= 0)
            error('%sPML{%d}.power_alpha must be a non-negative scalar.', direction, ii);
        end
    end

    % At this point, PML{ii} should contain 7 fields; check that no other fields exist
    if numel(fieldnames(PML{ii})) > 7
        field_names = setdiff(fieldnames(PML{ii}), {'npixels', 'power_sigma', 'sigma_max_over_omega', 'power_kappa', 'kappa_max', 'power_alpha', 'alpha_max_over_omega'});
        for jj = 1:numel(field_names)
            warning('%sPML{%d} contains unrecognized field ''%s''.', direction, ii, field_names{jj});
        end
    end
end

if PML{1}.npixels + PML{2}.npixels == 0
    PML = [];
end

end
