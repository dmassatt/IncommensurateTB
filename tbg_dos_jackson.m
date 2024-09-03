% Compute TBG density of states using KPM with Jackson smoothing. 
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
dE = 0.01;    % Energy grid spacing

load(['cheb_wgts_data/',filename]); % Load parameters and Chebyshev weights from file
cheb_wgts = cheb_wgts(1:p,:,:); % Truncate Chebyshev weights

E = (-E_range):dE:E_range; % Energy grid
nE = length(E);
[X,Y] = meshgrid(0:(N-1),0:(N-1));
X = X/N-1/2;
Y = Y/N-1/2;
X = X(:);
Y = Y(:);

[S,O] = meshgrid(1:2,1:2); % fix to just 1 if you want LDOS
S = S(:);
O = O(:);
ldos = zeros(size(E,2),size(X(:),1));

fprintf('begin loop\n')
for i = 1:size(X(:),1)
fprintf('%d / %d shift-loop\n',i,size(X(:),1))
    for j =1:4
        disp('Computing LDOS by Jackson KPM...')
        tic;
        Esc = E/(E_range+1);
        jackson_coeff = Cheb_JacksonCoeff(p-1);
        measure_weight = 1./sqrt(1 - Esc.^2);
        cheb_energy = Cheb_Eval(Esc, p-1);
        d = [.5 ones(1,p-1)];
        cheb_energy = diag(d)*cheb_energy;
        ldos(:,i) = ldos(:,i) + (((jackson_coeff.*cheb_wgts(:,i,j).') * cheb_energy) .*measure_weight)';
        disp(['Time=',num2str(toc)])
    end
end

dos = zeros(size(ldos(:,1)));
for i = 1:size(X(:),1)
    dos = dos + ldos(:,i);
end
dos = dos/(4*N^2); % normalize by discretization, # orbitals, # sheets

%addpath ~/Documents/MATLAB/export_fig/

figure(1);
plot(E, dos, '.-');
%export_fig('dos_jackson_1000.pdf');