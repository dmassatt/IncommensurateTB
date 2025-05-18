% Compute monolayer graphene density of states using KPM with Jackson smoothing. 
% Local Chebyshev weights <v|T_n(H)|v> computed using get_cheb_wgts_dos script.
%
% We load a 3-dim array cheb_wgts produced by that script. The first dimension
% the Chebyshev polynomial degree, the second orbital index.
% as the Hamiltonian is periodic, a single site is used to compute dos

% Input parameters
filename = 'r100_N4_p1000_dos.mat';
p = 501;        % Chebyshev degree
dE = 0.01;      % Energy grid spacing

L=100;
N = 2*L+1;

addpath('Tight_Binding')

E_range = 4;

E = (-E_range):dE:E_range; % Energy grid
nE = length(E);



dos = zeros(size(E,2),1);
j=1; % index

disp('Computing LDOS by Jackson KPM...')
tic;
    Esc = E/(E_range+1);
v = zeros(2*N^2,1);
v(1+(j-1)*N^2) = 1;
H = Generate_MonG(L,1); 
cheb_wgts=Cheb_LDoS_Weights(H, E_range, v, p-1);
jackson_coeff = Cheb_JacksonCoeff(p-1);
measure_weight = 1./sqrt(1 - Esc.^2);
cheb_energy = Cheb_Eval(Esc, p-1);
d = [.5 ones(1,p-1)];
cheb_energy = diag(d)*cheb_energy;
dos(:,j) = dos(:,j) + (((jackson_coeff.*cheb_wgts) * cheb_energy) .*measure_weight)';
disp(['Time=',num2str(toc)])


%addpath ~/Documents/MATLAB/export_fig/

figure(1);
plot(E, dos, '.-');
%export_fig('dos_jackson_1000.pdf');