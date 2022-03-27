%% Metalens: meta-atom
% Example of TiO2 hyperbolic metalens using SCSA-c
% 
% Use MESTI2S() to compute the transmission coefficient of meta-atom with different 
% ridge width and find meta-atoms satisfying 8 discrete ideal relative phase over 
% [0, 2pi)

%% System parameters

clear

% System parameters
n_air = 1; % Refractive index of air
n_silica = 1.46; % Refractive index of silica
n_TiO2 = 2.43; % Refractive index of TiO2
lambda = 532; % Free-space wavelength [nm]
dx = 13.3; % Grid size of system [nm]
w = 239.4; % Width of meta-atom cell [nm]
l = 600; % Thickness of meta-atom cell [nm]
ridge_hight = l; % Ridge height is the thickness of meta-atom cell.

%% General setup for  mesti2s()

% Setup input arguments for mesti2s(). 
syst.epsilon_L = n_silica^2; % Relative permittivity on the left hand side
syst.epsilon_R = n_air^2; % Relative permittivity on the right hand side
syst.wavelength = lambda; % Free-space wavelength [nm]
syst.dx = dx; % Grid size of system [nm]
syst.length_unit = 'nm'; %  Length unit
syst.yBC = 'periodic'; % Periodic boundary in y direction
in = {'left'}; % Specify input channel on the left.
out = {'right'}; % Specify output channel on the right.
opts.verbal = false; % Suppress output information.

% Plot refractive index profile of meta-atom
ridge_width = 79.8; % Ridge width of meta-atom [nm]
% Build permittivity for the meta-atom. 
% Please refer to the function build_epsilon_meta_atom.    
epsilon_meta_atom = build_epsilon_meta_atom(dx, n_air, n_TiO2, ridge_width, ridge_hight, w);
[ny, nx]= size(epsilon_meta_atom);
n_extra_for_plot = 10; % Extra pixels on side for plotting
% For ploting the space position
y = [0 ny*dx];
x = [-n_extra_for_plot*dx (nx+n_extra_for_plot)*dx]; 

figure
imagesc(x, y, [syst.epsilon_L*ones(ny,n_extra_for_plot), epsilon_meta_atom, 1*syst.epsilon_R*ones(ny,n_extra_for_plot)])
colormap(flipud(pink));
xlabel('{\itx} (nm)');
ylabel('{\ity} (nm)');
set(gca, 'fontsize', 15, 'FontName','Arial')
caxis([1 12])
text(670,130,'air','FontSize',20,'Rotation',90)
text(240,120,'TiO_2','FontSize',20)
text(-80,140,'silica','FontSize',20,'Rotation',90)
title(['Meta-atom'],'FontSize',20)

%% Transmission coefficient of meta-atom with different ridge width

ridge_width_list = 40:0.1:200; % List of ridge width: from 40 nm to 200 nm with 0.1 nm increment
t_list = zeros(1,size(ridge_width_list,2)); % Transmission coefficient list

% Looping over different ridge width
for ii =1:length(ridge_width_list)
    ridge_width = ridge_width_list(ii); % Ridge width of meta-atom [nm]
    syst.epsilon = build_epsilon_meta_atom(dx, n_air, n_TiO2, ridge_width, ridge_hight, w);
    % Call mesti2s() to calculate the scattering matrix.
    [S, stat] = mesti2s(syst, in, out, opts);        
    t_list(1,ii) = S(1,1);
end

phi0 = angle(t_list(ridge_width_list==194)); % Use the ridge width = 194 nm meta-atom as a phase reference.
rel_phi_over_pi_list = mod(angle(t_list)-phi0, 2*pi)/pi; % Relative phase over different ridge width

%%
% Plot the relative phase of meta-atom with different ridge width
figure
plot(ridge_width_list, rel_phi_over_pi_list, '-','linewidth', 2)
xlabel('Ridge width (nm)')
ylabel('Phase (\pi)')
xlim([40 200])
title('$\Phi - \Phi^0$', 'Interpreter','latex')
set(gca, 'fontsize', 20, 'FontName','Arial')
set(gca,'linewidth', 2)

%% Finding meta-atoms satisfying 8 discrete ideal relative phase over [0, 2pi)

ideal_rel_phase_over_pi_list = [linspace(0.25, 1.75, 7) 0]; % Have 8 discrete ideal relative phases over [0, 2pi).

% Find meta-atoms which are closeset to the ideal relative phase through nearest neighbor interpolation.
ind = interp1(rel_phi_over_pi_list,1:length(rel_phi_over_pi_list),ideal_rel_phase_over_pi_list,'nearest');
phi_over_pi_design_list = rel_phi_over_pi_list(ind); 
ridge_width_desgin_list = ridge_width_list(ind); 

% Print the relative phases and ridge widths.
fprintf(['Relative phase.(pi) %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f\n' ...
         'Ridge width....(nm) %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f\n'],...
         phi_over_pi_design_list,ridge_width_desgin_list);
% Save the phase list and the ridge width list.
save('meta_atom.mat','ridge_width_desgin_list','phi_over_pi_design_list')