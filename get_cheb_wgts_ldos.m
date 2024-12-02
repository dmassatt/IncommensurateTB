% Compute Chebyshev weights <v|T_n(H)|v> for use in local density of states
% calculuation

addpath("chebyshev");

p = 1000;           % # polynomials in Chebyshev expansion
E_range = 13;       % energy range
r_cut = 100;        % Radial truncation
theta = 6*pi/180;   % Rotation angle
filename = 'r100_p1000_ldos.mat';

X = -1/2;
Y = -1/2;

f = @(x) norm(x) < r_cut; % intralayer cut-off function

S = 1;
O = 1;

tstart = tic;

loc_s = graphene_init(theta,f,f,r_cut); % loc_s stores geometry, hopping function, system information.
v = loc_s.center_vector(S,O);
L = loc_s.sheet1.Lattice;

b = L*[X;Y];
disp('Generating Hamiltonian...')
tic;
%H = loc_s.MatrixShift(b);
H = GenerateH(theta,r_cut);
disp(['Time=',num2str(toc)])

disp('Generating Chebyshev weights...')
tic;
cheb_wgts = Cheb_LDoS_Weights(H, E_range, v, p-1).';
disp(['Time=',num2str(toc)])

disp(['Total time=',num2str(toc(tstart))])

save(['cheb_wgts_data/',filename], 'E_range', 'r_cut', 'theta', 'cheb_wgts');