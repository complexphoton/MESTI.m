%% Metalens: meta-atom
% Example of TiO2 hyperbolic metalens using SCSA-c
% 
% Use MESTI2S() to 
% 
% 1. Computing the transmission coefficient of meta-atom with different ridge 
% widths and find meta-atoms satisfying 8 discrete ideal relative phase over 
% [0, 2pi).
% 2. Scanning over ridge width and incident angle to get phase and amplitude 
% map of transmission coefficient.

%% System parameters
% Set up the parameters for the meta-atom system.

clear
addpath C:\Users\hcusc\Documents\USC_meeting\20220407\20220404

n_air    = 1;    % Refractive index of air
n_silica = 1.46; % Refractive index of silica
n_TiO2   = 2.43; % Refractive index of TiO2
lambda   = 532;  % Free-space wavelength [nm]
dx = lambda/40;  % Discretization grid size [nm]
w  = 18*dx;      % Width of meta-atom cell [nm]
l  = 600;        % Thickness of meta-atom cell [nm]
ridge_hight = l; % Ridge height is the thickness of meta-atom cell.

%% General setup for mesti2s()
% Set up input arguments for mesti2s(). 

syst.epsilon_L = n_silica^2; % Relative permittivity on the left hand side
syst.epsilon_R = n_air^2;    % Relative permittivity on the right hand side
syst.wavelength = lambda;    % Free-space wavelength [nm]
syst.dx = dx;                % Grid size of system [nm]
syst.length_unit = 'nm';     % Length unit
syst.yBC = 'periodic';       % Periodic boundary in y direction
in = {'left'};               % Specify input channel on the left.
out = {'right'};             % Specify output channel on the right.
opts.verbal = false;         % Suppress output information.

%% Structure of meta-atom
% Plot refractive index profile of a meta-atom with ridge with = 79.8 nm as 
% an illustration.

ridge_width = 79.8; % Ridge width of meta-atom [nm]

% Build permittivity for the meta-atom. 
% Please refer to the function build_epsilon_meta_atom.    
epsilon_meta_atom = build_epsilon_meta_atom(dx, n_air, n_TiO2, ridge_width, ridge_hight, w);
[ny, nx]= size(epsilon_meta_atom);

n_extra_for_plot = 10; % Extra pixels on side for plotting

% For plotting the space position
y = [0 ny*dx];
x = [-n_extra_for_plot*dx (nx+n_extra_for_plot)*dx]; 

figure
imagesc(x, y, [syst.epsilon_L*ones(ny,n_extra_for_plot), epsilon_meta_atom, 1*syst.epsilon_R*ones(ny,n_extra_for_plot)])
colormap(flipud(pink));
xlabel('{\itx} (nm)');
ylabel('{\ity} (nm)');
set(gca, 'fontsize', 15, 'FontName','Arial')
caxis([1 12])
xlim([x(1), x(2)])
ylim([y(1), y(2)])
text(670,130,'air','FontSize',20,'Rotation',90)
text(240,120,'TiO_2','FontSize',20)
text(-80,140,'silica','FontSize',20,'Rotation',90)
title(['Meta-atom'],'FontSize',20)

%% Transmission coefficient of meta-atom with different ridge widths
% In standard procedure of designing metasurface, people calculate the phase 
% of transmission coefficient of meta-atom with different parameters (for example, 
% ridge width). Then, use this information to design metasurface. Here following 
% the similar procedure, we calculate transmission coefficient of meta-atom looping 
% over different ridge widths.

ridge_width_list = 40:0.1:200; % List of ridge width: from 40 nm to 200 nm with 0.1 nm increment

t_list = zeros(1,size(ridge_width_list,2)); % Transmission coefficient list

% Looping over different ridge widths
for ii =1:length(ridge_width_list)
    ridge_width = ridge_width_list(ii); % Ridge width of meta-atom [nm]
    syst.epsilon = build_epsilon_meta_atom(dx, n_air, n_TiO2, ridge_width, ridge_hight, w);
    % Call mesti2s() to calculate the scattering matrix.
    [smatrix, stat] = mesti2s(syst, in, out, opts);

    % in = {'left'} and out = {'right'},
    % the smatrix = [t], where t is transmission coefficient from left to right.    
    t_list(1,ii) = smatrix(1,1);
end

phi0 = angle(t_list(ridge_width_list==194)); % Use the ridge width = 194 nm meta-atom as a phase reference.
rel_phi_over_pi_list = mod(angle(t_list)-phi0, 2*pi)/pi; % Relative phase over different ridge widths

% Plot the relative phase of meta-atom with different ridge widths
figure
plot(ridge_width_list, rel_phi_over_pi_list, '-','linewidth', 2)
xlabel('Ridge width (nm)')
ylabel('Phase (\pi)')
xlim([40 200])
title('$\Phi - \Phi^0$', 'Interpreter','latex')
set(gca, 'fontsize', 20, 'FontName','Arial')
set(gca,'linewidth', 2)

%% Finding meta-atoms satisfying 8 discrete ideal relative phases over [0, 2pi)
% In the design point of view, practically people would only use meta-atom with 
% discrete paratmeters to construct the metasurface. Typically, 8 different meta-atoms 
% covering the relative phases equally spacing over [0, 2pi) are used.

ideal_rel_phase_over_pi_list = [linspace(0.25, 1.75, 7) 0]; % Have 8 discrete equally spacing ideal relative phases over [0, 2pi).

% Find meta-atoms which are closest to the ideal relative phase through nearest neighbor interpolation.
ind = interp1(rel_phi_over_pi_list,1:length(rel_phi_over_pi_list),ideal_rel_phase_over_pi_list,'nearest');
phi_over_pi_design_list = rel_phi_over_pi_list(ind); 
ridge_width_design_list = ridge_width_list(ind); 

% Print the relative phases and ridge widths.
fprintf(['Relative phase.(pi) %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f\n' ...
         'Ridge width....(nm) %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f\n'],...
         phi_over_pi_design_list,ridge_width_design_list);
% Save the phase list and the ridge width list.
save('meta_atom.mat','ridge_width_design_list','phi_over_pi_design_list')

%% Phase and amplitude map of transmission coefficient of meta-atom with different ridge widths and incident angles
% In addition to normal incidence, the response of oblique incidence is also 
% important. With Bloch periodic boundary condition, the transmission coefficient 
% with different ridge widths and incident angles is calculated.

syst.yBC = 'Bloch'; % Bloch periodic boundary along transverse direction

ridge_width_list = 40:4:200; % List of ridge width: from 40 nm to 200 nm with 4 nm increment
theta_in_list = -89:1:89; % List of incident angle [degree]

k0dx = 2*pi/lambda*dx; % Dimensionless frequency k0dx

% Given theta, solve the corresponding kydx by the finite-difference dispersion
% Look at Eq. (S22) in the supplementary of the SCSA paper.
syms x
kydx_list = zeros(1,size(theta_in_list,2));
for jj = 1:size(theta_in_list,2)
    if theta_in_list(jj) == 0
        kydx_list(jj) = 0;    
    elseif theta_in_list(jj) > 0
        eqn = (k0dx)^2*n_silica^2 == 4*(sin(x/2))^2 + 4*(sin(x/tan(asin(1/n_silica*sind(theta_in_list(jj))))/2))^2;
        kydx_list(jj) = vpasolve(eqn,x,[0 2*pi]);
    elseif theta_in_list(jj) < 0        
        eqn = (k0dx)^2*n_silica^2 == 4*(sin(x/2))^2 + 4*(sin(x/tan(asin(1/n_silica*sind(theta_in_list(jj))))/2))^2;
        kydx_list(jj) = vpasolve(eqn,x,[-2*pi 0]);
    end
end

t_list = zeros(size(kydx_list,2), size(ridge_width_list,2));  % Transmission coefficient list 
% Row index for different incident and column index for different ridge widths

phi0_list = zeros(size(kydx_list,2), size(ridge_width_list,2)); % Reference phase list

% Looping over different ridge widths
for ii = 1:length(ridge_width_list)
    ridge_width = ridge_width_list(ii); % Ridge width of meta-atom [nm]
    syst.epsilon = build_epsilon_meta_atom(dx, n_air, n_TiO2, ridge_width, ridge_hight, w);

    % Looping over different incident angles
    for jj = 1:round(length(kydx_list))
        syst.ky_B = kydx_list(jj)/dx; % Bloch wave number in y direction
        % Call mesti2s() to calculate the scattering matrix.
        [smatrix, channels] = mesti2s(syst, in, out, opts); 

        % In some incident angles, there are more than one channel on the left.
        % Proper index should be chosen to extract the correct transmission coefficient for the phase/amplitude map.
        theta_inc = asind(sind(atand(channels.L.kydx_prop./channels.L.kxdx_prop))*n_silica); % Incident angle w.r.t. air of channels on the left
        ind = find(round(theta_inc)==theta_in_list(jj)); % Find the channel index whose incident angle users want.
        t_list(jj,ii) = smatrix(1,ind);

        % Choose the ridge width = 40 nm meta-atom as a phase reference.
        if ridge_width == 40
            phi0_list(jj,:) = angle(smatrix(1,ind)); 
        end
    end
end

% Plot the phase map of transmission coefficient over different ridge widths and incident angles.
figure
imagesc(ridge_width_list,theta_in_list, mod(angle(t_list)-phi0_list, 2*pi))
caxis([0, 2*pi]);
xlabel('Pillar width (nm)')
ylabel('\theta_{in} (degree)')
title('Phase')
colormap(twilight)
colorbar
hcb=colorbar; hcb.Ticks = [0 pi 2*pi]; hcb.TickLabels = {'0','\pi','2\pi'};
% Plot the amplitude map of transmission coefficient over different ridge widths and incident angles.

figure
imagesc(ridge_width_list,theta_in_list, abs(t_list))
caxis([0, 1]);
xlabel('Pillar width (nm)')
ylabel('\theta_{in} (degree)')
title('Amplitude')
colormap('hot')
colorbar