% Compute TBG local density of states using KPM with Jackson smoothing.
% Measure convergence with respect to p.
% Local Chebyshev weights <v|T_n(H)|v> computed using get_cheb_wgts_ldos script.
% We load a vector cheb_wgts of weights produced by that script.
%
% We generate a plot of the LDOS for different values of p, and compute
% the value of the LDOS at a specific point for different values of p.

% Input parameters
filename = 'r800_p4000_ldos.mat';
dE = 0.005;      % Energy grid spacing

ps = [500,1000,2000,4000];  % Polynomial degrees
E0 = -0.4;                  % Pick a specific energy to measure

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
    xlim([-2 1])

    measure_weight = 1./sqrt(1 - (E0/(E_range+1)).^2);
    cheb_energy = Cheb_Eval(E0/(E_range+1), p-1);
    d = [.5 ones(1,p-1)];
    cheb_energy = diag(d)*cheb_energy;
    ldos_val(i) = (((jackson_coeff.*cheb_wgts(1:p).') * cheb_energy) .*measure_weight)';
end
hold off


%export_fig('ldos_jackson_1000.pdf');