function [r_list_analytical,t_list_analytical] = fp_analytical(n_bg, n_slab, thickness, lambda_list)
%FP_ANALYTICAL calculate analytical reflection and transmission coefficient in 1D Fabry-Pérot (FP) etalon.
%
%   === Input Arguments ===
%   n_bg (numeric scalar, real or complex):
%       Refractive index of background material
%   n_slab (numeric scalar, real or complex):
%       Refractive index of dielectric slab
%   thickness (numeric scalar, real):
%       Thickness of the dielectric slab
%   lambda_list (numeric row vector):
%       List of wavelength
%   === Output Arguments ===
%   r_list_analytical (numeric row vector):
%       List of analytical reflection coefficient from left to left
%   t_list_analytical (numeric row vector):
%       List of analytical transmission coefficient from left to right

n_lambda = size(lambda_list,2); % Total number of wavelength to be calculated
r_list_analytical = zeros(1,n_lambda); % List of reflection coefficient
t_list_analytical = zeros(1,n_lambda); % List of transmission coefficient

% Looping over different wavelength to calculate reflection and transmission coefficient
for ii = 1:n_lambda    
    wavelength = lambda_list(ii); % wavelength 
    k0 = 2*pi/wavelength; % wavevector

    % According to the formula 7.1-63 in Fundamental of photonics by Saleh and Teich
    Mi = (1/(2*n_slab))*[n_bg+n_slab, n_slab-n_bg; n_slab-n_bg, n_bg+n_slab]; % Transfer matrix of the entrance boundary
    Mo = [exp(1i*n_slab*k0*thickness), 0; 0, exp(-1i*n_slab*k0*thickness)];  % Transfer matrix inside the slab
    Me = (1/(2*n_bg))*[n_slab+n_bg, n_bg-n_slab; n_bg-n_slab, n_slab+n_bg]; % Transfer matrix of the exit boundary
    M = Me*Mo*Mi; % Total transfer matrix
    % According to the formula 7.1-6 in Fundamental of photonics by Saleh and Teich, the relation between scattering matrix and transfer matrix 
    S = 1/M(2,2)*[M(1,1)*M(2,2)-M(1,2)*M(2,1), M(1,2) ; -M(2,1), 1];
    r_list_analytical(ii) = S(2,1); % Analytical reflection coefficient from left to left
    t_list_analytical(ii) = S(1,1); % Analytical transmission coefficient from left to right
end

end
