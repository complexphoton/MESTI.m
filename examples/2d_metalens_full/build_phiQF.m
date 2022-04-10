function [phiQF_L, phiQF_R] = build_phiQF(ny, N_L, N_R, w_t_over_dx_L, w_t_over_dx_R, use_Hann_Q)
%BUILD_PHIQF build the matrix phi*Q*F which is pre-processing step for SCSA-c.
%To learn more details, see Sec.5 in the supplementary of the SCSA paper.       
%
%   === Input Arguments ===
%   ny (numeric scalar, integer):
%       Number of pixels along the transverse direction
%   N_L (numeric scalar, integer):
%       Total number of channels on the left considering 
%   N_R (numeric scalar, integer):
%       Total number of channels on the right considering 
%   w_t_over_dx_L (numeric scalar, real):
%       Truncation window over grid size on the left
%   w_t_over_dx_R (numeric scalar, real):
%       Truncation window over grid size on the right
%   use_Hann_Q (logical scalar):
%       Whether to use Hann function as weight matrix Q.
%   === Output Arguments ===
%   phiQF_L (numeric matrix):
%       Matrix phi*Q*F on the left
%   phiQF_R (numeric matrix):
%       Matrix phi*Q*F on the right


% Truncation window in pixel
ny_tr_L = min([ny, 1 + round(w_t_over_dx_L)]);
ny_tr_R = min([ny, 1 + round(w_t_over_dx_R)]);

% Channel range to be considered
k_L = round((N_L-1)/2);
k_R = round((N_R-1)/2);

% Build the truncated version of phi_L*Q_L*F_L as a sparse matrix.
phiQF_L = spalloc(ny, N_L, ny_tr_L*N_L); % Pre-allocate the sparse matrix. 
for a_L = -k_L:1:k_L % Loop over input channels on the left.
    % Keep region along y direction (index m).
    m_center = ny + ny*a_L/N_L; % Peak location; not necessarily an integer
    m = round(m_center - ((ny_tr_L+1)/2)) + (1:ny_tr_L);  % In interval [1, 2*ny-1]
    m(m > ny) = m(m > ny) - ny; % Wrap around in y direction.

    if use_Hann_Q
        % Analytic expression for phi_L*Q_L*F_L
        % See Eq. (S36) and Eq. (S38) in the supplementary of the SCSA paper paper.       
        sin_m = sin((pi/ny*N_L)*m);
        term1 = (0.5*((-1)^a_L)/sqrt(ny*N_L))*sin_m./sin(pi*((m/ny)-(a_L/N_L)));
        term2 = (0.25*((-1)^(a_L-1))/sqrt(ny*N_L))*sin_m./sin(pi*((m/ny)-((a_L-1)/N_L)));
        term3 = (0.25*((-1)^(a_L+1))/sqrt(ny*N_L))*sin_m./sin(pi*((m/ny)-((a_L+1)/N_L)));
        % Handle situations where the denominator is zero.
        if int32(ny*a_L/N_L)*N_L == int32(ny*a_L); term1(m == (round(ny*a_L/N_L)+ny*(a_L<=0))) = 0.5*sqrt(N_L/ny); end
        if int32(ny*(a_L-1)/N_L)*N_L == int32(ny*(a_L-1)); term2(m == (round(ny*(a_L-1)/N_L)+ny*((a_L-1)<=0))) = 0.25*sqrt(N_L/ny); end
        if int32(ny*(a_L+1)/N_L)*N_L == int32(ny*(a_L+1)); term3(m == (round(ny*(a_L+1)/N_L)+ny*((a_L+1)<=0))) = 0.25*sqrt(N_L/ny); end
        phiQF_L(m, a_L+k_L+1) = term1 + term2 + term3;
    else
        % Analytic expression for phi_L*F_L
        % See Eq. (S36) in the supplementary of the SCSA paper paper.
        phiQF_L(m, a_L+k_L+1) = (((-1)^a_L)/sqrt(ny*N_L))*sin((pi/ny*N_L)*m)./sin(pi*((m/ny)-(a_L/N_L)));
        % Handle situations where the denominator is zero.        
        if int32(m_center)*N_L == int32(ny*(N_L+a_L)) 
            phiQF_L(round(m_center)-ny*(a_L>0), a_L+k_L+1) = sqrt(N_L/ny);
        end
    end
end

% Build the truncated version of phi_R*Q_R*F_R as a sparse matrix.
phiQF_R = spalloc(ny, N_R, ny_tr_R*N_R); % Pre-allocate the sparase matrix.
for a_R = -k_R:1:k_R  % Loop over input channels on the right.
    % Keep region along y direction (index m).
    m_center = ny + ny*a_R/N_R;  % Peak location; not necessarily an integer
    m = round(m_center - ((ny_tr_R+1)/2)) + (1:ny_tr_R); 
    m(m > ny) = m(m > ny) - ny; % Wrap around in y direction.

    if use_Hann_Q
        % Analytic expression for phi_R*Q_R*F_R
        % See Eq. (S36) and Eq. (S38) in the supplementary of the SCSA paper paper.       
        sin_m = sin((pi/ny*N_R)*m);
        term1 = (0.5*((-1)^a_R)/sqrt(ny*N_R))*sin_m./sin(pi*((m/ny)-(a_R/N_R)));
        term2 = (0.25*((-1)^(a_R-1))/sqrt(ny*N_R))*sin_m./sin(pi*((m/ny)-((a_R-1)/N_R)));
        term3 = (0.25*((-1)^(a_R+1))/sqrt(ny*N_R))*sin_m./sin(pi*((m/ny)-((a_R+1)/N_R)));
        % Handle situations where the denominator is zero.
        if int32(ny*a_R/N_R)*N_R == int32(ny*a_R); term1(m == (round(ny*a_R/N_R)+ny*(a_R<=0))) = 0.5*sqrt(N_R/ny); end
        if int32(ny*(a_R-1)/N_R)*N_R == int32(ny*(a_R-1)); term2(m == (round(ny*(a_R-1)/N_R)+ny*((a_R-1)<=0))) = 0.25*sqrt(N_R/ny); end
        if int32(ny*(a_R+1)/N_R)*N_R == int32(ny*(a_R+1)); term3(m == (round(ny*(a_R+1)/N_R)+ny*((a_R+1)<=0))) = 0.25*sqrt(N_R/ny); end
        phiQF_R(m, a_R+k_R+1) = term1 + term2 + term3;
    else
        % Analytic expression for phi_R*F_R        
        % See Eq. (S36) in the supplementary of the SCSA paper paper.
        phiQF_R(m, a_R+k_R+1) = (((-1)^a_R)/sqrt(ny*N_R))*sin((pi/ny*N_R)*m)./sin(pi*((m/ny)-(a_R/N_R)));
        % Handle situations where the denominator is zero.                
        if int32(m_center)*N_R == int32(ny*(N_R+a_R))
            phiQF_R(round(m_center)-ny*(a_R>0), a_R+k_R+1) = sqrt(N_R/ny);
        end
    end
end

end
