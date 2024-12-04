% Compute TBG local density of states using KPM with Jackson smoothing.
% Measure convergence with respect to eta.
% Local Chebyshev weights <v|T_n(H)|v> computed using get_cheb_wgts_ldos script.
% We load a vector cheb_wgts of weights produced by that script.
%
% We generate a plot of the LDOS for different values of eta, and compute
% the value of the LDOS at a specific point for different values of eta.

% Input parameters
filename = 'r800_p4000_ldos.mat';
p = 4000;     % Chebyshev degree
m = 6;        % Order of method with respect to broadening parameter eta
dE = 0.005;    % Energy grid spacing

etas = 0.02./2.^(0:3);   % Broadening parameters
E0 = -0.4; % Pick specific energy to measure

addpath('hodc','hodc/kernels');

load(['cheb_wgts_data/',filename]); % Load parameters and Chebyshev weights from file

E = (-E_range):dE:E_range; % Energy grid

ldos_val = zeros(size(etas));
figure(3);
% Compute local densities of states
for i=1:length(etas)
    eta = etas(i);
    ldos = hodc_ldos(m, eta, p, E/E_range, cheb_wgts(1:p));

    plot(E, ldos, '.-'); hold on
    xlim([-2, 1])

    ldos_val(i) = hodc_ldos(m, eta, p, E0/E_range, cheb_wgts(1:p));
end
hold off

legend(string(num2cell(etas)));
xlabel('E')
ylabel('LDOS(E)')
set(gca,'fontsize',20)


%export_fig('dos_hodc_1000.pdf');