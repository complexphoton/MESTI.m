function [r_list_analytical, t_list_analytical] = dbr_analytical(n_bg, n_1, n_2, d_1, d_2, N_pair, lambda_list)
%DBR_ANALYTICAL calculate analytical reflection and transmission coefficient in 1D distributed Bragg reflector (DBR).
%
%   === Input Arguments ===
%   n_bg (numeric scalar, real or complex):
%       Refractive index of background material
%   n_1 (numeric scalar, real or complex):
%       Refractive index of material 1
%   n_2 (numeric scalar, real or complex):
%       Refractive index of material 2
%   d_1 (numeric scalar, real):
%       Thickness of the material 1 in a pair
%   d_2 (numeric scalar, real):
%       Thickness of the material 2 in a pair
%   N_pair (numeric scalar, integer):
%       Number of pair in DBR
%   lambda_list (numeric row vector):
%       List of wavelength
%   === Output Arguments ===
%   r_list_analytical (numeric row vector):
%       List of analytical reflection coefficient from left to left
%   t_list_analytical (numeric row vector):
%       List of analytical transmission coefficient from left to right

n_lambda = size(lambda_list,2); % Number of angular frequencies to be calculated
r_list_analytical = zeros(1,n_lambda); % List of reflection coefficient
t_list_analytical = zeros(1,n_lambda); % List of transmission coefficient

for ii = 1:n_lambda    
    lambda = lambda_list(ii); % wavelength 
    k0 = 2*pi/lambda; % wavevector

    % According to the formula 7.1-63 in Fundamental of photonics by Saleh and Teich
    Mi = (1/(2*n_1))*[(n_1+n_bg), (n_1-n_bg); (n_1-n_bg), (n_1+n_bg)]; % Transfer matrix of the entrance boundary
    Mo = (1/(4*n_1*n_2))*...
         [(n_1+n_2)*exp(1i*n_2*k0*d_2), (n_1-n_2)*exp(-1i*n_2*k0*d_2); ...
         (n_1-n_2)*exp(1i*n_2*k0*d_2), (n_1+n_2)*exp(-1i*n_2*k0*d_2)]*...
         [(n_2+n_1)*exp(1i*n_1*k0*d_1), (n_2-n_1)*exp(-1i*n_1*k0*d_1); ...
         (n_2-n_1)*exp(1i*n_1*k0*d_1), (n_2+n_1)*exp(-1i*n_1*k0*d_1)]; % Transfer matrix of in one pair of DBR  
    Me = (1/(4*n_bg*n_2))*...
         [(n_bg+n_2)*exp(1i*n_2*k0*d_2), (n_bg-n_2)*exp(-1i*n_2*k0*d_2);...
         (n_bg-n_2)*exp(1i*n_2*k0*d_2),(n_bg+n_2)*exp(-1i*n_2*k0*d_2)]*...
         [(n_2+n_1)*exp(1i*n_1*k0*d_1), (n_2-n_1)*exp(-1i*n_1*k0*d_1); ...
         (n_2-n_1)*exp(1i*n_1*k0*d_1), (n_2+n_1)*exp(-1i*n_1*k0*d_1)]; % Transfer matrix of Nth segment with a boundary into the exit boundary
    M = Me*(Mo)^(N_pair-1)*Mi; % Total transfer matrix
    % According to the formula 7.1-6 in Fundamental of photonics by Saleh and Teich, the relation between scattering matrix and transfer matrix     
    S = 1/M(2,2)*[M(1,1)*M(2,2)-M(1,2)*M(2,1), M(1,2) ; -M(2,1), 1];
    r_list_analytical(ii) = S(2,1); % Analytical reflection coefficient from left to left
    t_list_analytical(ii) = S(1,1); % Analytical transmission coefficient from left to right
end

end
