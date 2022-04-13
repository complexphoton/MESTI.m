function epsilon_meta_atom = build_epsilon_meta_atom(dx, n_bg, n_ridge, ridge_width, ridge_height, w)

%BUILD_EPSILON_META_ATOM generate relative permittivity profile of meta-atom with subpixel smoothing.
%
%   === Input Arguments ===
%   dx (numeric scalar, real):
%       Grid size               
%   n_bg (numeric scalar, real or complex):
%       Refractive index of background material
%   n_ridge (numeric scalar, real or complex):
%       Refractive index of ridge material
%   ridge_width (numeric scalar, real):
%       Width of the ridge
%   ridge_height (numeric scalar, real):
%       Height of the ridge
%   w (numeric scalar, real):
%       Width of meta-atom cell
%   === Output Arguments ===
%   epsilon_meta_atom (numeric matrix, real or complex):
%       Discretized relative permittivity profile of meta-atom

% Number of pixel for the meta-atom 
nx = ceil(ridge_height/dx);    
ny = ceil(w/dx);

ny_ridge = ridge_width/dx; % Number of pixels to be assigned as ridge

% Setting the permittivity profile with subpixel smoothing
epsilon_meta_atom = n_bg^2*ones(ny, nx);
if mod(ny,2) == 0 % When ny is even number
    pixel_ridge_ratio_y = mod(ny_ridge,2)/2;
    pixel_bg_ratio_y = 1-pixel_ridge_ratio_y;
    epsilon_meta_atom((ny/2-floor(ny_ridge/2)+1):(ny/2+floor(ny_ridge/2)),:) = (n_ridge)^2;
    epsilon_meta_atom([(ny/2-floor(ny_ridge/2)+1-1), (ny/2+floor(ny_ridge/2)+1)],:) = (n_ridge)^2*pixel_ridge_ratio_y+(n_bg)^2*pixel_bg_ratio_y;
    
    pixel_bg_ratio_x = nx - ridge_height/dx;
    pixel_ridge_ratio_x = 1 - pixel_bg_ratio_x;
    epsilon_meta_atom((ny/2-floor(ny_ridge/2)+1):(ny/2+floor(ny_ridge/2)),end) = (n_ridge)^2*pixel_ridge_ratio_x+(n_bg)^2*pixel_bg_ratio_x;
    epsilon_meta_atom([(ny/2-floor(ny_ridge/2)+1-1),(ny/2+floor(ny_ridge/2)+1)],end) = (n_ridge)^2*pixel_ridge_ratio_x*pixel_ridge_ratio_y+(n_bg)^2*(1-pixel_ridge_ratio_x*pixel_ridge_ratio_y);
elseif ny_ridge > 1 % When ny is odd number and ny_ridge > 1
    pixel_ridge_ratio_y = mod(ny_ridge-1,2)/2;
    pixel_bg_ratio_y = 1-pixel_ridge_ratio_y;
    epsilon_meta_atom(((ny+1)/2-floor((ny_ridge-1)/2)):((ny+1)/2+floor((ny_ridge-1)/2)),:) = (n_ridge)^2;
    epsilon_meta_atom([((ny+1)/2-floor((ny_ridge-1)/2)-1), ((ny+1)/2+floor((ny_ridge-1)/2)+1)],:) = (n_ridge)^2*pixel_ridge_ratio_y+(n_bg)^2*pixel_bg_ratio_y;
    
    pixel_bg_ratio_x = nx - ridge_height/dx;
    pixel_ridge_ratio_x = 1 - pixel_bg_ratio_x;
    epsilon_meta_atom(((ny+1)/2-floor((ny_ridge-1)/2)):((ny+1)/2+floor((ny_ridge-1)/2)),end) = (n_ridge)^2*pixel_ridge_ratio_x+(n_bg)^2*pixel_bg_ratio_x;
    epsilon_meta_atom([((ny+1)/2-floor((ny_ridge-1)/2)-1), ((ny+1)/2+floor((ny_ridge-1)/2)+1)],end) = (n_ridge)^2*pixel_ridge_ratio_x*pixel_ridge_ratio_y+(n_bg)^2*(1-pixel_ridge_ratio_x*pixel_ridge_ratio_y);
else % When ny is odd number and ny_ridge <= 1
    pixel_ridge_ratio_y = ny_ridge;
    pixel_bg_ratio_y = 1-pixel_ridge_ratio_y;
    epsilon_meta_atom((ny+1)/2,:) = (n_ridge)^2*pixel_ridge_ratio_y+(n_bg)^2*pixel_bg_ratio_y;

    pixel_bg_ratio_x = nx - ridge_height/dx;
    pixel_ridge_ratio_x = 1 - pixel_bg_ratio_x;
    epsilon_meta_atom((ny+1)/2,end) = (n_ridge)^2*pixel_ridge_ratio_x*pixel_ridge_ratio_y+(n_bg)^2*(1-pixel_ridge_ratio_x*pixel_ridge_ratio_y);
end