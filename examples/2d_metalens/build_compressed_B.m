function phiQF = build_compressed_B(ny, M, ny_window, use_Hann_window)
%BUILD_PHIQF build the compressed matrix phi*Q*F for APF-c.
%   See supplementary section 5 of the APF paper.
%
%   === Input Arguments ===
%   ny (positive integer scalar):
%       Number of pixels along the transverse direction.
%   M (positive integer scalar):
%       Number of channels. Must be no larger than ny.
%   ny_window (positive integer scalar):
%       Number of pixels within the truncation window. Must be between 1 and ny.
%   use_Hann_window (logical scalar):
%       Whether to use the Hann window function as the weight Q.
%
%   === Output Arguments ===
%   phiQF (ny-by-M sparse matrix):
%       Matrix phi*Q*F truncated to ny_window pixels around the peaks.

if M > ny
    error('M should not be larger than ny.');
end

if ny_window > ny
    warning('ny_window = %d is larger than ny = %d; reducing to ny.');
    ny_window = ny;
end

% Range of channel index a
a_min = round((M-1)/2);

% Build the truncated version of phi*Q*F as a sparse matrix.
m_L_window = 1:ny_window;
phiQF = spalloc(ny, M, ny_window*M); % Pre-allocate the sparse matrix. 
for a = -a_min:1:a_min % Loop over input channels on the left.
    % Find the indices m for pixels within the truncation window
    m_center = ny + ny*a/M; % Peak location; not necessarily an integer
    m = round(m_center - ((ny_window+1)/2)) + m_L_window; % In interval [1, 2*ny-1]
    m(m > ny) = m(m > ny) - ny; % Wrap around in y direction.

    % Build phi*Q*F within this window
    if use_Hann_window
        % See Eq. (S36) and Eq. (S38) in the supplementary of the APF paper.
        sin_m = sin((pi/ny*M)*m);
        term1 = (0.5*((-1)^a)/sqrt(ny*M))*sin_m./sin(pi*((m/ny)-(a/M)));
        term2 = (0.25*((-1)^(a-1))/sqrt(ny*M))*sin_m./sin(pi*((m/ny)-((a-1)/M)));
        term3 = (0.25*((-1)^(a+1))/sqrt(ny*M))*sin_m./sin(pi*((m/ny)-((a+1)/M)));
        % Handle situations where the denominator is zero.
        if int32(ny*a/M)*M == int32(ny*a); term1(m == (round(ny*a/M)+ny*(a<=0))) = 0.5*sqrt(M/ny); end
        if int32(ny*(a-1)/M)*M == int32(ny*(a-1)); term2(m == (round(ny*(a-1)/M)+ny*((a-1)<=0))) = 0.25*sqrt(M/ny); end
        if int32(ny*(a+1)/M)*M == int32(ny*(a+1)); term3(m == (round(ny*(a+1)/M)+ny*((a+1)<=0))) = 0.25*sqrt(M/ny); end
        phiQF(m, a+a_min+1) = term1 + term2 + term3;
    else
        % See Eq. (S36) in the supplementary of the APF paper.
        phiQF(m, a+a_min+1) = (((-1)^a)/sqrt(ny*M))*sin((pi/ny*M)*m)./sin(pi*((m/ny)-(a/M)));
        % Handle situations where the denominator is zero.
        if int32(m_center)*M == int32(ny*(M+a))
            phiQF(round(m_center)-ny*(a>0), a+a_min+1) = sqrt(M/ny);
        end
    end
end

end
