% Compute TBG local density of states using KPM with Jackson smoothing.
% Measure convergence with respect to eta.
% Local Chebyshev weights <v|T_n(H)|v> computed using get_cheb_wgts_ldos script.
% We load a vector cheb_wgts of weights produced by that script.
%
% We generate a plot of the LDOS for different values of eta, and compute
% the value of the LDOS at a specific point for different values of eta.

addpath("chebyshev");
addpath("hodc");
% Input parameters
filename = 'r1600_p8000_ldos.mat';
p = 8000;     % Chebyshev degree
m = 6;        % Order of method with respect to broadening parameter eta
dE = 0.005;    % Energy grid spacing

etas = 0.08 ./2.^(0:7);   % Broadening parameters

% Value just above 0.01
new_value = 0.012;
new_value2 =  0.011;
nv1 = 0.009;
nv2 = 0.008;
nv3 = 0.007;
nv4 = 0.006;
% Add the new value to the array and sort it
etas = sort([etas, new_value,new_value2,nv1,nv2,nv3,nv4]);
E0 = -0.43; % Pick specific energy to measure

addpath('hodc','hodc/kernels');

load(['cheb_wgts_data/',filename]); % Load parameters and Chebyshev weights from file

%E = (-E_range):dE:E_range; % Energy grid

ldos_val = zeros(size(etas));
%figure(3);
% Compute local densities of states
for i=1:length(etas)
    eta = etas(i);
    %ldos = hodc_ldos(m, eta, p, E/E_range, cheb_wgts(1:p));

    %plot(E, ldos, '.-'); hold on
    %xlim([-1.5, 0.5])

    ldos_val(i) = hodc_ldos(m, eta, p, E0/E_range, cheb_wgts(1:p));
end

% Combine eta and ldos_val into a two-column matrix
data = [etas(:), ldos_val(:)];

% Save the data to a .dat file
filename = 'E0_ldos_values_hodc_r1600_p8000.dat';
dlmwrite(filename, data, 'delimiter', '\t', 'precision', '%.6f');


hold off

%legend(string(num2cell(etas)));
%xlabel('E')
%ylabel('LDOS(E)')
%set(gca,'fontsize',20)

%exportgraphics(gcf, 'ldos_hodc_8000.png', 'Resolution', 300)
%export_fig('dos_hodc_1000.pdf');