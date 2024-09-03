% Compute TBG density of states using high-order delta Chebyshev method. 
% Local Chebyshev weights <v|T_n(H)|v> computed using get_cheb_wgts script.
%
% We load a 3-dim array cheb_wgts produced by that script. The first dimension
% the Chebyshev polynomial degree, the second indexes shifts, and the third
% indexes sheet/orbital combinations.

% Input parameters (must be same as in dos_compute.m)
filename = 'r100_N4_p1000.mat';
p = 501;      % Chebyshev degree
% N = 4;        % # shifts per dimension
% E_range = 13; % Energy range covering entire spectrum
m = 6;        % Order of method with respect to broadening parameter eta
eta = 0.03;   % Broadening parameter
dE = 0.01;    % Energy grid spacing

load(['cheb_wgts_data/',filename]); % Load parameters and Chebyshev weights from file
cheb_wgts = cheb_wgts(1:p,:,:); % Truncate Chebyshev weights

E = (-E_range):dE:E_range; % Energy grid
nE = length(E);

% Compute local densities of states correspondning to shifts and sheet/orbital
% combinations
ldos = hodc_ldos(m, eta, p, E/E_range, reshape(cheb_wgts, p, N^2 * 4));
ldos = reshape(ldos, nE, N^2, 4);
ldos = sum(ldos, 3); % Sum over sheet/orbital combinations

% Compute density of states
dos = sum(ldos, 2)/(4*N^2); % Average over shifts

%addpath ~/Documents/MATLAB/export_fig/

figure(1);
plot(E, dos, '.-');
%export_fig('dos_hodc_1000.pdf');