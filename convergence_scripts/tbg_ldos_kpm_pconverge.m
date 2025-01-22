% Compute TBG local density of states using KPM with Jackson smoothing.
% Measure convergence with respect to p.
% Local Chebyshev weights <v|T_n(H)|v> computed using get_cheb_wgts_ldos script.
% We load a vector cheb_wgts of weights produced by that script.
%
% We generate a plot of the LDOS for different values of p, and compute
% the value of the LDOS at a specific point for different values of p.

% Input parameters
filename = 'r800_p16000_ldos.mat';
dE = 0.005;      % Energy grid spacing

ps = [500,1000,2000,4000];  % Polynomial degrees

addpath('../kpm')
load(['../cheb_wgts_data/',filename]); % Load parameters and Chebyshev weights from file

E = (-E_range):dE:E_range; % Energy grid

figure(5);
% Compute local densities of states
for i=1:length(ps)
    p = ps(i);
    ldos = kpm_ldos(p, E/E_range, cheb_wgts(1:p));

    plot(E, ldos, '.-'); hold on
    xlim([-2 1])
end
hold off

legend(string(num2cell(ps)));
xlabel('E')
ylabel('LDOS(E)')
set(gca,'fontsize',15)