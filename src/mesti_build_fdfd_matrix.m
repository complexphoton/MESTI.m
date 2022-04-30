function [A, is_symmetric_A, xPML, yPML] = mesti_build_fdfd_matrix(eps_or_inv_eps, k0dx, xBC, yBC, xPML, yPML, use_UPML)
%MESTI_BUILD_FDFD_MATRIX The finite-difference frequency-domain operator in 2D.
%   A = mesti_build_fdfd_matrix(epsilon, k0dx, xBC, yBC, xPML, yPML, use_UPML)
%   returns A as a sparse matrix representing wave operator
%      [- (d/dx)^2 - (d/dy)^2 - k0^2*epsilon(x,y)]*(dx^2)
%   for the Ez component of transverse-magnetic (TM) waves (Hx, Hy, Ez),
%   discretized on a square grid with grid size dx through center difference.
%
%   A = mesti_build_fdfd_matrix(inv_epsilon, k0dx, xBC, yBC, xPML, yPML, use_UPML)
%   returns A as a sparse matrix representing wave operator
%      [- (d/dx)*(inv_epsilon{2})*(d/dx) - (d/dy)*(inv_epsilon{1})*(d/dy) - k0^2]*(dx^2)
%   for the Hz component of transverse-electric (TE) waves (Ex, Ey, Hz).
%
%   === Input Arguments ===
%   eps_or_inv_eps (matrix or two-element cell array; required):
%      When eps_or_inv_eps is a matrix, TM polarization is considered, and
%      epsilon = eps_or_inv_eps is a ny-by-nx numeric matrix (real or complex)
%      discretizing the relative permittivity profile on integer sites.
%         When eps_or_inv_eps is two-element cell array, TE polarization is
%      considered, and inv_epsilon = eps_or_inv_eps with
%         inv_epsilon{1}: ny_d-by-nx matrix discretizing inv(epsilon(x,y))_xx
%         inv_epsilon{2}: ny-by-nx_d matrix discretizing inv(epsilon(x,y))_yy
%      where inv(epsilon(x,y)) is a diagonal tensor.
%         Here, ny_d is the number of grid points of Ex ~ dHz/dy on half-integer
%      sites in y, and nx_d is the number of grid points of Ey ~ dHz/dx on half-
%      integer sites in x; they depend on the boundary condition: n_d = n for
%      Bloch periodic, DirichletNeumann, and NeumannDirichlet boundary
%      conditions; n_d = n + 1 for Dirichlet boundary condition, and n_d = n - 1
%      for Neumann boundary condition.
%   k0dx (numeric scalar, real or complex; required):
%      Normalized frequency k0*dx = (2*pi/vacuum_wavelength)*dx.
%   xBC (character vector or numeric scalar; required):
%      Boundary condition (BC) at the two ends in x direction, effectively
%      specifying Ez(m,n) or Hz(m,n) at n=0 and n=nx+1 which are one pixel
%      beyond our computation domain. For character inputs, the options are:
%         'periodic'         - f(m,n+nx) = f(m,n)
%         'Dirichlet'        - f(m,0) = f(m,nx+1) = 0
%         'Neumann'          - f(m,0) = f(m,1); f(m,nx+1) = f(m,nx)
%         'DirichletNeumann' - f(m,0) = 0; f(m,nx+1) = f(m,nx)
%         'NeumannDirichlet' - f(m,0) = f(m,1); f(m,nx+1) = 0
%      where f = Ez or Hz depending on the polarization.
%         One can also specify 'PEC', 'PMC', 'PECPMC', or 'PMCPEC'. For TM
%      polarization, PEC is equivalent to Dirichlet for which Ez = 0 at the
%      boundary, and PMC is equivalent to Neumann for which dEz/dx or dEz/dy = 0
%      at the boundary. For TE polarization, PMC is equivalent to Dirichlet for
%      which Hz = 0 at the boundary, and PEC is equivalent to Neumann for which
%      dHz/dx or dHz/dy = 0 at the boundary.
%         When xBC is a numeric scalar, the Bloch periodic boundary condition is
%      used with f(m,n+nx) = f(m,n)*exp(1i*xBC); in other words, xBC =
%      kx_B*nx*dx = kx_B*p where kx_B is the Bloch wave number and p = nx*dx is
%      the periodicity in x.
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
%            Note this is within syst.epsilon or syst.inv_epsilon.
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
%      (ny*nx)-by-(ny*nx) sparse matrix representing the 2D FDFD operator.
%   is_symmetric_A (logical scalar):
%      Whether matrix A is symmetric or not.
%   xPML (two-element cell array or []):
%      PML parameters used on the low and high sides of x direction, if any.
%   yPML (two-element cell array or []):
%      PML parameters used on the low and high sides of y direction, if any.

% Check input arguments
% xBC and yBC are checked in build_ddx_1d(); xPML and yPML are checked in set_PML_params()
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
    if numel(inv_epsilon) ~= 2
        error('Input argument inv_epsilon must be a two-element cell array.');
    elseif ~(ismatrix(inv_epsilon{1}) && isnumeric(inv_epsilon{1}))
        error('inv_epsilon{1} must be a numeric matrix.');
    elseif ~(ismatrix(inv_epsilon{2}) && isnumeric(inv_epsilon{2}))
        error('inv_epsilon{2} must be a numeric matrix.');
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

% Number of grid points in y and x
if use_TM
    [ny, nx] = size(epsilon);
else
    [ny_d_eps, nx] = size(inv_epsilon{1}); % inv_epsilon_xx
    [ny, nx_d_eps] = size(inv_epsilon{2}); % inv_epsilon_yy
end

% Nothing to do if the size is zero
nxy = nx*ny;
if nxy == 0
    A = sparse(0, 0);
    return
end

% Compute background permittivity for PML
if use_TM
    epsilon_bg_x = real([mean(epsilon(:,1)), mean(epsilon(:,end))]);
    epsilon_bg_y = real([mean(epsilon(1,:)), mean(epsilon(end,:))]);
else
    % Note that nx_d and ny_d can be zero with Neumann BC, in which case we can't take the mean.
    % If nx_d = 0 or ny_d = 0, there must be no PML in that direction, so epsilon_bg is not needed.
    if nx_d_eps ~= 0
        epsilon_bg_x = real(1./[mean(inv_epsilon{2}(:,1)), mean(inv_epsilon{2}(:,end))]);
    else
        epsilon_bg_x = zeros(2,1);
    end
    if ny_d_eps ~= 0
        epsilon_bg_y = real(1./[mean(inv_epsilon{1}(1,:)), mean(inv_epsilon{1}(end,:))]);
    else
        epsilon_bg_y = zeros(2,1);
    end
end

% Set default values for PML parameters
xPML = set_PML_params(xPML, k0dx, epsilon_bg_x, 'x');
yPML = set_PML_params(yPML, k0dx, epsilon_bg_y, 'y');

% Convert PEC/PMC to Dirichlet/Neumann based on whether TE or TM polarization is used
xBC = convert_BC(xBC, use_TM);
yBC = convert_BC(yBC, use_TM);

% Build the first derivatives in 1D, such that ddx*f = df/dx, ddy*f = df/dy (ignoring dx factor)
% The next derivative is given by -ddx'. For example, d^2f/dx^2 = (-ddx')*ddx*f.
% Note that sx and sx_d are column vectors.
% Later we will use Kronecker outer product to go from 1D to 2D.
[ddx, sx, sx_d] = build_ddx_1d(nx, xBC, xPML, 'x');
[ddy, sy, sy_d] = build_ddx_1d(ny, yBC, yPML, 'y');

% Number of grid points for df/dx and df/fy
nx_d = size(ddx, 1);
ny_d = size(ddy, 1);

if ~use_TM
    % Check if size(inv_epsilon{1}) and size(inv_epsilon{2}) are compatible with the boundary conditions
    if nx_d_eps ~= nx_d
        error('nx_d = size(inv_epsilon{2},2) = %d but should be %d instead given nx = size(inv_epsilon{1},2) = %d and the boundary condition in x.', nx_d_eps, nx_d, nx);
    elseif ny_d_eps ~= ny_d
        error('ny_d = size(inv_epsilon{1},1) = %d but should be %d instead given ny = size(inv_epsilon{2},1) = %d and the boundary condition in y.', ny_d_eps, ny_d, ny);
    end
    inv_epsilon_xx = spdiags(inv_epsilon{1}(:), 0, nx*ny_d, nx*ny_d);
    inv_epsilon_yy = spdiags(inv_epsilon{2}(:), 0, nx_d*ny, nx_d*ny);
end

% Convert to diagonal matrices
inv_sx_d = spdiags(1./sx_d, 0, nx_d, nx_d);
inv_sy_d = spdiags(1./sy_d, 0, ny_d, ny_d);

% Build the operators
% TM: A = [- (d/dx)^2 - (d/dy)^2 - k0^2*epsilon(x,y)]*(dx^2)
% TE: A = [- (d/dx)*(inv_epsilon_yy)*(d/dx) - (d/dy)*(inv_epsilon_xx)*(d/dy) - k0^2]*(dx^2)
% Note we cancel the minus sign of (-ddx') with the overall minus sign.
if ~use_UPML
    % Stretched-coordinate PML (SC-PML); here, 1/s and 1/s_d are both multiplied onto the gradient.
    inv_sx = spdiags(1./sx, 0, nx, nx);
    inv_sy = spdiags(1./sy, 0, ny, ny);
    if use_TM
        A = kron(inv_sx*(ddx')*inv_sx_d*ddx, speye(ny)) ...
          + kron(speye(nx), inv_sy*(ddy')*inv_sy_d*ddy) ...
          - spdiags((k0dx^2)*epsilon(:), 0, nxy, nxy);
    else
        A = kron(inv_sx*(ddx'), speye(ny)) * inv_epsilon_yy * kron(inv_sx_d*ddx, speye(ny)) ...
          + kron(speye(nx), inv_sy*(ddy')) * inv_epsilon_xx * kron(speye(nx), inv_sy_d*ddy) ...
          - (k0dx^2)*speye(nxy);
    end
else
    % Uniaxial PML (UPML); here, A_UPML = (sx*sy)*A_SCPML
    s_xy = sy.*reshape(sx, 1, nx); % use implicit expansion
    sx = spdiags(sx, 0, nx, nx);
    sy = spdiags(sy, 0, ny, ny);
    if use_TM
        A = kron((ddx')*inv_sx_d*ddx, sy) ...
          + kron(sx, (ddy')*inv_sy_d*ddy) ...
          - spdiags((k0dx^2)*(s_xy(:).*epsilon(:)), 0, nxy, nxy);
    else
        A = kron(ddx', sy) * inv_epsilon_yy * kron(inv_sx_d*ddx, speye(ny)) ...
          + kron(sx, ddy') * inv_epsilon_xx * kron(speye(nx), inv_sy_d*ddy) ...
          - spdiags((k0dx^2)*s_xy(:), 0, nxy, nxy);
    end
end

if (isnumeric(xBC) && xBC ~= 0 && xBC ~= pi) || (isnumeric(yBC) && yBC ~= 0 && yBC ~= pi)
    % Bloch periodic boundary condition with ka != 0 or pi breaks the symmetry of A
    is_symmetric_A = false;
elseif (isempty(xPML) && isempty(yPML)) || use_UPML
    is_symmetric_A = true;
else
    % SC-PML also breaks the symmetry of A
    is_symmetric_A = false;
end

end


function BC_new = convert_BC(BC, use_TM)
% Convert PEC/PMC to Dirichlet/Neumann based on whether TE or TM polarization is used
BC_new = BC;
if use_TM
    % For TM waves, we consider Ez(x,y), so PEC => Dirichlet, PMC => Neumann
    if strcmpi(BC, 'PEC')
        BC_new = 'Dirichlet';
    elseif strcmpi(BC, 'PMC')
        BC_new = 'Neumann';
    elseif strcmpi(BC, 'PECPMC')
        BC_new = 'DirichletNeumann';
    elseif strcmpi(BC, 'PMCPEC')
        BC_new = 'NeumannDirichlet';
    end
else
    % For TE waves, we consider Hz(x,y), so PMC => Dirichlet, PEC => Neumann
    if strcmpi(BC, 'PMC')
        BC_new = 'Dirichlet';
    elseif strcmpi(BC, 'PEC')
        BC_new = 'Neumann';
    elseif strcmpi(BC, 'PMCPEC')
        BC_new = 'DirichletNeumann';
    elseif strcmpi(BC, 'PECPMC')
        BC_new = 'NeumannDirichlet';
    end
end

end


function [grad, s, s_d] = build_ddx_1d(n, BC, PML, direction)
% grad: first-derivative matrix
% s: coordinate-stretching factor for f and d^2f/dx^2 (on integer sites)
% s_d: coordinate-stretching factor for df/dx (on half-integer sites)

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

% Build the first-derivative matrix
% f = [f(1), ..., f(n)].' on integer sites
% df = (df/dx)*dx on half-integer sites
if strcmpi(BC, 'Bloch')
    % f(n+1) = f(1)*exp(1i*ka); f(0) = f(n)*exp(-1i*ka)
    % grad*f = df = [df(1.5), ..., df(n+0.5)].'
    grad = spdiags([ones(n,1),-ones(n,1),exp(1i*ka)*ones(n,1)], [1,0,1-n], n, n);
elseif strcmpi(BC, 'Dirichlet') % Dirichlet on both sides
    % f(0) = f(n+1) = 0
    % grad*f = df = [df(0.5), ..., df(n+0.5)].'
    grad = spdiags([ones(n,1),-ones(n,1)], [0,-1], n+1, n);
elseif strcmpi(BC, 'Neumann') % Neumann on both sides
    % df(0.5) = df(n+0.5) = 0
    % grad*f = df = [df(1.5), ..., df(n-0.5)].'
    grad = spdiags([ones(n-1,1),-ones(n-1,1)], [1,0], n-1, n);
elseif strcmpi(BC, 'DirichletNeumann') % Dirichlet on the low side, Neumann on the high side
    % f(0) = 0; df(n+0.5) = 0
    % grad*f = df = [df(0.5), ..., df(n-0.5)].'
    grad = spdiags([ones(n,1),-ones(n,1)], [0,-1], n, n);
elseif strcmpi(BC, 'NeumannDirichlet') % Neumann on the low side, Dirichlet on the high side
    % df(0.5) = 0; f(n+1) = 0
    % grad*f = df = [df(1.5), ..., df(n+0.5)].'
    grad = spdiags([ones(n,1),-ones(n,1)], [1,0], n, n);
else
    error('Input argument %sBC = %s is not a supported option.', direction, BC);
end

n_d = size(grad,1); % number of sites for df/dx

% Coordinate-stretching factor: s(x) = kappa(x) + sigma(x)/(alpha(u) - i*omega)
s = ones(n,1); % s-factor for f and d^2f/dx^2 (on integer sites)
s_d = ones(n_d,1); % s-factor for df/dx (on half-integer sites)

% No coordinate stretching if no PML is specified
if isempty(PML)
    return
end

% Number of PML pixels on the low and high sides
npixels = [PML{1}.npixels, PML{2}.npixels];

% Cannot have more PML pixels than the number of pixels
% In our implementation, the conductivity goes to zero at one pixel before PML. So there needs to be at least one more pixel in addition to PML.
if sum(npixels) >= n
    error('Total number of PML pixels = %d + %d = %d should be less than the number of pixels = %d in %s direction.', npixels(1), npixels(2), sum(npixels), n, direction);
end

% Below, u(x) is a function that goes linearly from u(x)=0 at one site before PML to u(x)=1 at the "end of PML".
% Note that where the "end of PML" is depends on the boundary condition.
% Let index i=0 be one site before PML, so u(i=0)=0.
% Then, i=1 is the first site of PML and i=npixels is the last site of PML we explicitly simulate. But we do not set u(i=npixels)=1.
% For Dirichlet BC, the BC is such that f=0 at i=npixels+1, so we let the end of PML be i=npixels+1, with an effective PML thickness of npixels+1 pixels.
% For Neumann BC, the BC is such that df=0 at i=npixels+0.5, so we let the end of PML be i=npixels+0.5, with an effective PML thickness of npixels+0.5 pixels.
% For periodic BC, the BC is such that df(n+0.5)=df(0.5)*exp(1i*ka), and we let the end of PML be i=npixels+0.5 so that if PML is used on both sides, we will have u(0.5)=u(n+0.5)=1 so the conductivity profile is continuous across the edge if sigma_max is the same on both sides.

% Construct u(x) and their corresponding indices for df/dx on half-integer sites
if strcmpi(BC, 'Bloch')
    npixels_effective = [npixels(1)+0.5, npixels(2)+0.5];
    if npixels(1)*npixels(2) == 0
        warning('Bloch periodic boundary condition is used with a single-sided PML in %s direction; transmission through PML will only undergo single-pass attenuation.', direction);
    end
    % no s-factor for df(0.5) since we only consider df(n+0.5)
    u_PML_d = {((1:npixels(1))-0.5)/npixels_effective(1), ((1:(npixels(2)+1))-0.5)/npixels_effective(2)};
    ind_PML_d = {fliplr(1:npixels(1)), (n+1)-fliplr(1:(npixels(2)+1))};
elseif strcmpi(BC, 'Dirichlet') % Dirichlet on both sides
    npixels_effective = [npixels(1)+1, npixels(2)+1];
    u_PML_d = {((1:(npixels(1)+1))-0.5)/npixels_effective(1), ((1:(npixels(2)+1))-0.5)/npixels_effective(2)};
    ind_PML_d = {fliplr(1:(npixels(1)+1)), (n+2)-fliplr(1:(npixels(2)+1))};
elseif strcmpi(BC, 'Neumann') % Neumann on both sides
    npixels_effective = [npixels(1)+0.5, npixels(2)+0.5];
    u_PML_d = {((1:npixels(1))-0.5)/npixels_effective(1), ((1:npixels(2))-0.5)/npixels_effective(2)};
    ind_PML_d = {fliplr(1:npixels(1)), n-fliplr(1:npixels(2))};
elseif strcmpi(BC, 'DirichletNeumann') % Dirichlet on the low side, Neumann on the high side
    npixels_effective = [npixels(1)+1, npixels(2)+0.5];
    u_PML_d = {((1:(npixels(1)+1))-0.5)/npixels_effective(1), ((1:npixels(2))-0.5)/npixels_effective(2)};
    ind_PML_d = {fliplr(1:(npixels(1)+1)), (n+1)-fliplr(1:npixels(2))};
elseif strcmpi(BC, 'NeumannDirichlet') % Neumann on the low side, Dirichlet on the high side
    npixels_effective = [npixels(1)+0.5, npixels(2)+1];
    u_PML_d = {((1:npixels(1))-0.5)/npixels_effective(1), ((1:(npixels(2)+1))-0.5)/npixels_effective(2)};
    ind_PML_d = {fliplr(1:npixels(1)), (n+1)-fliplr(1:(npixels(2)+1))};
else
    error('Input argument %sBC = ''%s'' is not a supported option.', direction, BC);
end

% Construct u(x) and their corresponding indices for f and d^2f/dx^2 on integer sites
u_PML = {(1:npixels(1))/npixels_effective(1), (1:npixels(2))/npixels_effective(2)};
ind_PML = {fliplr(1:npixels(1)), (n+1)-fliplr(1:npixels(2))};

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
        % Note s and s_d should be sampled directly from the polynomial functions; if we use the polynomial for one and linear interpolation for the other, the performance will be much worse.
        s(ind_PML{ii}) = func_s(u_PML{ii}).';  % column vector
        s_d(ind_PML_d{ii}) = func_s(u_PML_d{ii}).';  % column vector
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
