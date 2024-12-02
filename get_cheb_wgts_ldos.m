% Compute Chebyshev weights <v|T_n(H)|v> for use in local density of states
% calculuation

addpath("chebyshev");

p = 1000;           % # polynomials in Chebyshev expansion
E_range = 13;       % energy range
r_cut = 100;        % Radial truncation
theta = 6*pi/180;   % Rotation angle
filename = 'r100_p1000_ldos.mat';

S = 1; % sheet number for LDoS
O = 1; % orbital index for LDoS

tstart = tic;

disp('Generating Hamiltonian...')
tic;

[H,sheet_orbital] = GenerateH(theta,r_cut);
v = zeros(size(H,1),1);
v(sheet_orbital(S,O)) = 1;
disp(['Time=',num2str(toc)])

disp('Generating Chebyshev weights...')
tic;
cheb_wgts = Cheb_LDoS_Weights(H, E_range, v, p-1).';
disp(['Time=',num2str(toc)])

disp(['Total time=',num2str(toc(tstart))])

save(['cheb_wgts_data/',filename], 'E_range', 'r_cut', 'theta', 'cheb_wgts');