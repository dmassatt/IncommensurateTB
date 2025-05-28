% Compute TBG local density of states using high-order delta-Chebyshev method, for
% decreasing values of eta, and plot self-convergence error.
% Local Chebyshev weights <v|T_n(H)|v> computed using get_cheb_wgts_ldos script.
% We load a vector cheb_wgts of weights produced by that script.

% Input parameters
filename = 'r800_p16000_ldos.mat';
p = 8000;     % Chebyshev degree
m = 6;        % Order of method with respect to broadening parameter eta
dE = 0.005;    % Energy grid spacing

etas = 0.02./1.5.^(0:5);   % Broadening parameters

addpath('../hodc','../hodc/kernels');

load(['../cheb_wgts_data/',filename]); % Load parameters and Chebyshev weights from file

E = (-E_range):dE:E_range; % Energy grid

figure(3);
% Compute local densities of states
for i=1:length(etas)
    eta = etas(i);
    ldos = hodc_ldos(m, eta, p, E/E_range, cheb_wgts(1:p));

    plot(E, ldos, '.-'); hold on
    xlim([-2, 1])
end
hold off

legend(string(num2cell(etas)));
xlabel('E')
ylabel('LDOS(E)')
set(gca,'fontsize',20)