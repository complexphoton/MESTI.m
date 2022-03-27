function [epsilon, x0_list, y0_list, r0_list] = build_epsilon_disorder(W, L, r_min, r_max, min_sep, nummber_density, rng_seed, dx, epsilon_scat, epsilon_bg)
%   === Input Arguments ===
%   W       = width of the scattering region
%   L       = thickness of the scattering region
%   r_min   = minimal radius of the cylinders
%   r_max   = maximal radius of the cylinders
%   min_sep = minimal separation between cylinders
%   nummber_density = number density of the cylinders
%   rng_seed = random number generator seed
%   dx       = discretization grid size
%   epsilon_scat = relative permittivity of the cylinders
%   epsilon_bg   = relative permittivity of the background

N_scatterers = round(nummber_density*L*W);
nx = round(L/dx);
ny = round(W/dx);

% Edges of the box
x1 = 0;
x2 = L;
y1 = 0;
y2 = W;

% Randomly pick the radii and locations of N_scatterers non-overlapping cylinders, and then generate the discretized permittivity profile
rng(rng_seed); % set random number generator seed
x0_list = zeros(N_scatterers, 1);
y0_list = zeros(N_scatterers, 1);
r0_list = zeros(N_scatterers, 1);
epsilon = epsilon_bg*ones(ny,nx);
for ii = 1:N_scatterers
    % Radius of the cylinder
    r0 = r_min + rand*(r_max-r_min);

    % min/max values for the coordinates of the cylinder center
    % We forbid the cylinder to go over the y=0 or y=W boundary; this is to simplify the code below (so we don't need to handle the periodicity in y)
    x_min = x1 + r0;
    x_max = x2 - r0;
    y_min = y1 + r0;
    y_max = y2 - r0;

    % Keep randomly picking the cylinder center until we get one that doesn't overlap with the previous cylinders
    ind_prev = 1:(ii-1);
    found = false;
    while ~found
        x0 = x_min + rand*(x_max-x_min);
        y0 = y_min + rand*(y_max-y_min);
        if min((x0_list(ind_prev)-x0).^2+(y0_list(ind_prev)-y0).^2-(r0_list(ind_prev)+(r0+min_sep)).^2) <= 0
            found = false;
        else
            found = true;
        end
    end
    x0_list(ii) = x0;
    y0_list(ii) = y0;
    r0_list(ii) = r0;

    % Identify indices of the interiors of the cylinder within a local box surrounding (x0,y0) with edge length = 2*r0
    % Note that we don't do subpixel smoothing here for simplicity
    n0 = 0.5 + x0/dx; % location of (x0,y0) in terms of index (n from 1 to nx, m from 1 to ny)
    m0 = 0.5 + y0/dx;
    dn = r0/dx;
    n1 = max([1,  round(n0 - dn)]);
    n2 = min([nx, round(n0 + dn)]);
    m1 = max([1,  round(m0 - dn)]);
    m2 = min([ny, round(m0 + dn)]);
    x_local = ((n1-0.5):n2)*dx;
    y_local = ((m1-0.5):m2)*dx;
    [X_local, Y_local] = meshgrid(x_local,y_local);
    [m_local, n_local] = find((((X_local-x0).^2 + (Y_local-y0).^2)) <= (r0^2));

    % Set the relative permittivity of the cylinder
    ind = sub2ind([ny, nx], m_local+(m1-1), n_local+(n1-1)); % convert to linear indices in the full matrix
    epsilon(ind) = epsilon_scat;
end

end
