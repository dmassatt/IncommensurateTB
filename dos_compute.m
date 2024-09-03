% This first block of code computes the DoS

% assume existence of loc_s
r_cut = 100; % radial truncation of system

f = @(x) norm(x) < r_cut; % intralayer cut-off function
theta = 6*pi/180;

addpath("chebyshev");

P = 500; % chebyshev order

N = 4; % discretization of shifts

E_range = 13; % energy range
E = (-E_range):.01:E_range; % needs to encompass entire spectrum
[X,Y] = meshgrid(0:(N-1),0:(N-1));
X = X/N-1/2;
Y = Y/N-1/2;
X = X(:);
Y = Y(:);

[S,O] = meshgrid(1:2,1:2); % fix to just 1 if you want LDOS
S = S(:);
O = O(:);
ldos = zeros(size(E,2),size(X(:),1));
cheb_wgts = zeros(P+1,size(X(:),1),4);

loc_s = graphene_init(theta,f,f,r_cut); % loc_s stores geometry, hopping function, system information.
fprintf('begin loop\n')
for i = 1:size(X(:),1)
fprintf('%d / %d shift-loop\n',i,size(X(:),1))
    for j =1:4

        v = loc_s.center_vector(S(j),O(j));
        if S(j) == 1
            L = loc_s.sheet1.Lattice;
        else
            L = loc_s.sheet2.Lattice;
        end

        b = L*[X(i);Y(i)];
        disp('Generating Hamiltonian for this shift...')
        tic;
        H = loc_s.MatrixShift(b);
        disp(['Time=',num2str(toc)])

        disp('Computing LDOS by KPM...')
        tic;
        [ldostmp,cheb_wgts(:,i,j)] = Cheb_LDoS(H, E_range, v, P, E/(E_range+1));
        ldos(:,i) = ldos(:,i) + ldostmp';
        disp(['Time=',num2str(toc)])
    end
end

dos = zeros(size(ldos(:,1)));
for i = 1:size(X(:),1)
    dos = dos + ldos(:,i);
end
dos = dos/(4*N^2); % normalize by discretization, # orbitals, # sheets

%addpath ~/Documents/MATLAB/export_fig/

%% display DoS
figure(1);
plot(E,dos);
Range = 1; % Energy is displayed on [E_F-Range, E_F+Range]
E_F = -.6; % rough position of Fermi energy
%axis([-Range+E_F,Range+E_F,0,max(dos(abs(E-E_F)<Range))])
%export_fig('dos2.png');

%% Display all LDoS on same plot

figure(2);
plot(E,ldos);
Range = 1;
E_F = -.6; % rough position of Fermi energy
ldos_cut= ldos(abs(E-E_F)<Range,:);
%axis([-Range+E_F,Range+E_F,0, max(ldos_cut(:)) ])
%export_fig('ldos2.png');
