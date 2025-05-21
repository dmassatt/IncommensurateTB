% Compute graphene local density of states using high-order delta-Chebyshev method, for
% decreasing values of eta.
% Local Chebyshev weights <v|T_n(H)|v> computed using get_cheb_wgts_ldos script.
% We load a vector cheb_wgts of weights produced by that script.

% Input parameters
filename = 'graphene_L800_p16000_ldos.mat';
p = 4000;     % Chebyshev degree
m = 6;        % Order of method with respect to broadening parameter eta
dE = 0.002;    % Energy grid spacing

etas = 0.2./1.5.^(0:5);   % Broadening parameters

addpath('../hodc','../hodc/kernels');

load(['cheb_wgts_data/',filename]); % Load parameters and Chebyshev weights from file

E = (-E_range):dE:E_range; % Energy grid

figure(3);
% Compute local densities of states
for i=1:length(etas)
    eta = etas(i);
    % ldos = hodc_ldos(m, eta, p, E/E_range, cheb_wgts(1:p));
    ldos = hodc_ldos(m, eta, p, E, E_range, cheb_wgts(1:p));

    plot(E, ldos, '.-'); hold on
    xlim([-4, 4])
end
ldos_true = zeros(length(E),1);
for i=1:length(E)
    ldos_true(i) = graphene_analytic(E(i),1);
end
plot(E,ldos_true,'k--','LineWidth',2); % Analytic result
hold off

legend(string(num2cell(etas)));
xlabel('E')
ylabel('LDOS(E)')
set(gca,'fontsize',20)