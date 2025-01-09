% Compute TBG local density of states using KPM with Jackson smoothing.
% Measure convergence with respect to p.
% Local Chebyshev weights <v|T_n(H)|v> computed using get_cheb_wgts_ldos script.
% We load a vector cheb_wgts of weights produced by that script.
%
% We generate a plot of the LDOS for different values of p, and compute
% the value of the LDOS at a specific point for different values of p.

addpath("chebyshev");

% Input parameters
filename = 'r1600_p8000_ldos.mat';
dE = 0.005;      % Energy grid spacing

ps = [500,1000,2000,3000,4000,5000,6000,7000,8000];  % Polynomial degrees
E0 = -0.43;                  % Pick a specific energy to measure

load(['cheb_wgts_data/',filename]); % Load parameters and Chebyshev weights from file

E = (-E_range):dE:E_range; % Energy grid

disp('Computing LDOS by Jackson KPM...')

Esc = E/(E_range+1);

ldos_val = zeros(size(ps));
figure(4);
for i=1:length(ps)
    p = ps(i);
    jackson_coeff = Cheb_JacksonCoeff(p-1);
    measure_weight = 1./sqrt(1 - Esc.^2);
    cheb_energy = Cheb_Eval(Esc, p-1);
    d = [.5 ones(1,p-1)];
    cheb_energy = diag(d)*cheb_energy;
    ldos = (((jackson_coeff.*cheb_wgts(1:p).') * cheb_energy) .*measure_weight)';

    plot(E, ldos, '.-'); hold on
    xlim([-1.5 0.5])

    measure_weight = 1./sqrt(1 - (E0/(E_range+1)).^2);
    cheb_energy = Cheb_Eval(E0/(E_range+1), p-1);
    d = [.5 ones(1,p-1)];
    cheb_energy = diag(d)*cheb_energy;
    ldos_val(i) = (((jackson_coeff.*cheb_wgts(1:p).') * cheb_energy) .*measure_weight)';
end


% Combine eta and ldos_val into a two-column matrix
data = [ps(:), ldos_val(:)];

% Save the data to a .dat file
filename = 'E0_ldos_values_kpm_r1600.dat';
dlmwrite(filename, data, 'delimiter', '\t', 'precision', '%.6f');

hold off
%exportgraphics(gcf, 'ldos_kpm_8000.png', 'Resolution', 300)
%export_fig('ldos_jackson_8000.pdf');