% Compute TBG local density of states using KPM with Jackson smoothing. 
% Local Chebyshev weights <v|T_n(H)|v> computed using get_cheb_wgts_ldos script.
% We load a vector cheb_wgts of weights produced by that script.

% Input parameters
filename = 'r800_p4000_ldos.mat';
p = 4000;       % Chebyshev degree
dE = 0.005;      % Energy grid spacing

load(['cheb_wgts_data/',filename]); % Load parameters and Chebyshev weights from file

E = (-E_range):dE:E_range; % Energy grid

disp('Computing LDOS by Jackson KPM...')
tic;
Esc = E/(E_range+1);
jackson_coeff = Cheb_JacksonCoeff(p-1);
measure_weight = 1./sqrt(1 - Esc.^2);
cheb_energy = Cheb_Eval(Esc, p-1);
d = [.5 ones(1,p-1)];
cheb_energy = diag(d)*cheb_energy;
ldos = (((jackson_coeff.*cheb_wgts(1:p).') * cheb_energy) .*measure_weight)';
disp(['Time=',num2str(toc)])

%addpath ~/Documents/MATLAB/export_fig/

figure(1);
plot(E, ldos, '.-');
xlim([-2 1])
%export_fig('ldos_jackson_1000.pdf');