function [epsilon, inv_epsilon, x0_list, y0_list, r0_list] = build_epsilon_disorder(W, L, r_min, r_max, min_sep, number_density, rng_seed, dx, epsilon_scat, epsilon_bg, build_TM, build_TE, yBC)
%   Generate a random collection of cylinders and build the relative permittivity profile.
%   === Input Arguments ===
%   W       = width of the scattering region
%   L       = thickness of the scattering region
%   r_min   = minimal radius of the cylinders
%   r_max   = maximal radius of the cylinders
%   min_sep = minimal separation between cylinders
%   number_density = number density of the cylinders
%   rng_seed = random number generator seed
%   dx       = discretization grid size
%   epsilon_scat = relative permittivity of the cylinders
%   epsilon_bg   = relative permittivity of the background
%   build_TM = whether to build epsilon for TM polarization
%   build_TE = whether to build 1/epsilon for TE polarization
%   yBC = boundary condition in y

N_scatterers = round(number_density*L*W);
nx_Ez = round(L/dx);
ny_Ez = round(W/dx);

if build_TM
    epsilon = epsilon_bg*ones(ny_Ez,nx_Ez);
else
    epsilon = [];
end
if build_TE
    nx_Hz = nx_Ez + 1;
    %    For periodic, Bloch periodic, and PMCPEC boundary conditions in y,
    % ny_Hz = ny_Ez, and all of the sites on y_{n+0.5} are half a pixel
    % after the corresponding sites on y_n.
    %    For PEC boundary condition in y, ny_Hz = ny_Ez + 1, and the sites
    % on y_{n+0.5} start from half a pixel before the first site of y_n and
    % end on half a pixel after the last site of y_n.
    %    For PMC boundary condition in y, ny_Hz = ny_Ez - 1, and the sites
    % on y_{n+0.5} start from half a pixel after the first site of y_n and
    % end on half a pixel before the last site of y_n.
    %    For PECPMC boundary condition in y, ny_Hz = ny_Ez, and all of the
    % sites on y_{n+0.5} are half a pixel before the corresponding sites on
    % y_n.
    if ismember(lower(yBC), lower({'Bloch', 'periodic', 'PMCPEC'}))
        ny_Hz = ny_Ez;
        dm_Hz = 0;
    elseif strcmpi(yBC, 'PEC')
        ny_Hz = ny_Ez + 1;
        dm_Hz = 1;
    elseif strcmpi(yBC, 'PMC')
        ny_Hz = ny_Ez - 1;
        dm_Hz = 0;
    elseif strcmpi(yBC, 'PECPMC')
        ny_Hz = ny_Ez;
        dm_Hz = 1;
    else
        error('yBC = ''%s'' is not a supported option.', yBC);
    end
    inv_epsilon = cell(1,2);
    inv_epsilon{1} = (1/epsilon_bg)*ones(ny_Ez,nx_Hz); % (1/epsilon(x,y))_xx
    inv_epsilon{2} = (1/epsilon_bg)*ones(ny_Hz,nx_Ez); % (1/epsilon(x,y))_yy
else
    inv_epsilon = [];
end

% Edges of the box
x1 = 0;
x2 = L;
y1 = 0;
y2 = W;

if N_scatterers > 3000
    % Store the indices of the existing cylinders by their approximate locations, to facilitate overlap checking when N_scatterers is large
    use_cell = true;
    cell_size = 2*r_max + min_sep;
    nx_cell = ceil(L/cell_size);
    ny_cell = ceil(W/cell_size);
    ind_scatterers = cell(ny_cell, nx_cell);
else
    use_cell = false;
end

% Randomly pick the radii and locations of N_scatterers non-overlapping cylinders, and then generate the discretized permittivity profile
rng(rng_seed); % set random number generator seed
x0_list = zeros(N_scatterers, 1);
y0_list = zeros(N_scatterers, 1);
r0_list = zeros(N_scatterers, 1);
for ii = 1:N_scatterers
    % Radius of the cylinder
    r0 = r_min + rand*(r_max-r_min);

    if use_cell
        % Distance for overlap checking
        dist_check = r0 + min_sep + r_max;
    end

    % min/max values for the coordinates of the cylinder center
    % We forbid the cylinder to go over the y=0 or y=W boundary; this is to simplify the code below (so we don't need to handle the periodicity in y)
    x_min = x1 + r0;
    x_max = x2 - r0;
    y_min = y1 + r0;
    y_max = y2 - r0;

    % Keep randomly picking the cylinder center until we get one that doesn't overlap with the existing cylinders
    found = false;
    while ~found
        % Center of the new cylinder
        x0 = x_min + rand*(x_max-x_min);
        y0 = y_min + rand*(y_max-y_min);

        if use_cell
            % Indices of cells we need to check for overlap
            n1 = max([1      , round(0.5 + (x0 - dist_check)/cell_size)]);
            n2 = min([nx_cell, round(0.5 + (x0 + dist_check)/cell_size)]);
            m1 = max([1      , round(0.5 + (y0 - dist_check)/cell_size)]);
            m2 = min([ny_cell, round(0.5 + (y0 + dist_check)/cell_size)]);

            % Indices of cylinders we need to check for overlap
            ind_check = cell2mat(reshape(ind_scatterers(m1:m2,n1:n2), 1, []));

            % Check those cylinders for overlap
            if isempty(ind_check) || min((x0_list(ind_check)-x0).^2+(y0_list(ind_check)-y0).^2-(r0_list(ind_check)+(r0+min_sep)).^2) > 0
                found = true;
            end
        else
            % Check all existing cylinders for overlap
            if ii==1 || min((x0_list(1:(ii-1))-x0).^2+(y0_list(1:(ii-1))-y0).^2-(r0_list(1:(ii-1))+(r0+min_sep)).^2) > 0
                found = true;
            end
        end
    end
    x0_list(ii) = x0;
    y0_list(ii) = y0;
    r0_list(ii) = r0;

    if use_cell
        % Store the index of this cylinder by its approximate location
        n0 = min([nx_cell, max([1, round(0.5 + x0/cell_size)])]);
        m0 = min([ny_cell, max([1, round(0.5 + y0/cell_size)])]);
        ind_scatterers{m0, n0} = [ind_scatterers{m0, n0}, ii];
    end

    % Set the relative permittivity of the cylinder
    dn = r0/dx;
    if build_TM
        % Identify indices of the interiors of the cylinder within a local box surrounding (x0,y0) with edge length = 2*r0
        % Note that we don't do subpixel smoothing here for simplicity
        % Ez and epsilon are located at (x,y) = (n-0.5, m-0.5)*dx with n from 1 to nx_Ez, m from 1 to ny_Ez
        n0 = 0.5 + x0/dx; % location of (x0,y0) in terms of index
        m0 = 0.5 + y0/dx;
        n1 = max([1,  round(n0 - dn)]);
        n2 = min([nx_Ez, round(n0 + dn)]);
        m1 = max([1,  round(m0 - dn)]);
        m2 = min([ny_Ez, round(m0 + dn)]);
        x_local = ((n1-0.5):n2)*dx;
        y_local = ((m1-0.5):m2)*dx;
        [X_local, Y_local] = meshgrid(x_local,y_local);
        [m_local, n_local] = find((((X_local-x0).^2 + (Y_local-y0).^2)) <= (r0^2));
        ind = sub2ind([ny_Ez, nx_Ez], m_local+(m1-1), n_local+(n1-1)); % convert to linear indices in the full matrix
        epsilon(ind) = epsilon_scat;
    end
    if build_TE
        % inv_epsilon{1} = (1/epsilon)_xx is located at (x,y) = (n-1, m-0.5)*dx with n from 1 to nx_Hz = nx_Ez + 1, m from 1 to ny_Ez
        n0 =  1  + x0/dx; % location of (x0,y0) in terms of index
        m0 = 0.5 + y0/dx;
        n1 = max([1,  round(n0 - dn)]);
        n2 = min([nx_Hz, round(n0 + dn)]);
        m1 = max([1,  round(m0 - dn)]);
        m2 = min([ny_Ez, round(m0 + dn)]);
        x_local = ((n1:n2)-1)*dx;
        y_local = ((m1-0.5):m2)*dx;
        [X_local, Y_local] = meshgrid(x_local,y_local);
        [m_local, n_local] = find((((X_local-x0).^2 + (Y_local-y0).^2)) <= (r0^2));
        ind = sub2ind([ny_Ez, nx_Hz], m_local+(m1-1), n_local+(n1-1)); % convert to linear indices in the full matrix
        inv_epsilon{1}(ind) = 1/epsilon_scat;
        % Note that inv_epsilon{1}(:,1) and inv_epsilon{1}(:,end) are half outside the scattering region, so they need to be set separately given epsilon_L and epsilon_R.

        % inv_epsilon{2} = (1/epsilon)_yy is located at x = (n-0.5)*dx with n from 1 to nx_Ez, and
        %   y = m*dx for periodic, Bloch periodic, PMC, and PMCPEC boundary conditions in y
        %   y = (m-1)*dx for PEC and PECPMC boundary conditions in y
        % with m from 1 to ny_Hz
        n0 = 0.5 + x0/dx; % location of (x0,y0) in terms of index
        m0 = dm_Hz + y0/dx;
        n1 = max([1,  round(n0 - dn)]);
        n2 = min([nx_Ez, round(n0 + dn)]);
        m1 = max([1,  round(m0 - dn)]);
        m2 = min([ny_Hz, round(m0 + dn)]);
        x_local = ((n1-0.5):n2)*dx;
        y_local = ((m1:m2)-dm_Hz)*dx;
        [X_local, Y_local] = meshgrid(x_local,y_local);
        [m_local, n_local] = find((((X_local-x0).^2 + (Y_local-y0).^2)) <= (r0^2));
        ind = sub2ind([ny_Hz, nx_Ez], m_local+(m1-1), n_local+(n1-1)); % convert to linear indices in the full matrix
        inv_epsilon{2}(ind) = 1/epsilon_scat;
    end
end

end
