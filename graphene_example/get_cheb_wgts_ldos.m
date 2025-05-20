% Compute Chebyshev weights <v|T_n(H)|v> for monolayer graphene local density of
% states calculation

addpath("../kpm")
addpath("../Tight_Binding");
% maxNumCompThreads(1);

p = 16000;      % # polynomials in Chebyshev expansion
E_range = 3.5;  % energy range
L = 800;       % Radial truncation
filename = 'graphene_L800_p16000_ldos.mat';

tstart = tic;

fprintf('Generating LDOS Cheb wgts for L = %d, p = %d\n', L, p);
tic;

H = Generate_MonG(L,1);
v = zeros(size(H,1),1);
v(1) = 1;
fprintf('Time: %f s\n', toc);
fprintf('Total memory used: %d MB\n\n', sum([whos().bytes])/1e6);

disp('Generating Chebyshev weights...')
tic;
cheb_wgts = Cheb_LDoS_Weights(H, E_range, v, p-1).';

fprintf('Time: %f s\n', toc);
fprintf('Total memory used: %d MB\n\n', sum([whos().bytes])/1e6);

save(['cheb_wgts_data/',filename], 'E_range', 'L', 'cheb_wgts');
