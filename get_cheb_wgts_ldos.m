% Compute Chebyshev weights <v|T_n(H)|v> for use in local density of states
% calculuation

addpath("chebyshev");
maxNumCompThreads(1);

p = 16000;          % # polynomials in Chebyshev expansion
E_range = 13;       % energy range
r_cut = 3200;        % Radial truncation
theta = 6*pi/180;   % Rotation angle
filename = 'r3200_p16000_ldos.mat';

S = 1; % sheet number for LDoS
O = 1; % orbital index for LDoS

tstart = tic;

fprintf('Generating LDOS Cheb wgts for r = %d, p = %d, theta = %f\n', r_cut, p, theta);
fprintf('Initial memory used: %d MB\n\n', sum([whos().bytes])/1e6);
disp('Generating Hamiltonian...')
tic;

[H,sheet_orbital] = GenerateH(theta,r_cut);
v = zeros(size(H,1),1);
v(sheet_orbital(S,O)) = 1;
fprintf('Time: %f s\n', toc);
fprintf('Total memory used: %d MB\n\n', sum([whos().bytes])/1e6);

disp('Generating Chebyshev weights...')
tic;
cheb_wgts = Cheb_LDoS_Weights(H, E_range, v, p-1).';

fprintf('Time: %f s\n', toc);
fprintf('Total memory used: %d MB\n\n', sum([whos().bytes])/1e6);

save(['cheb_wgts_data/',filename], 'E_range', 'r_cut', 'theta', 'cheb_wgts');
