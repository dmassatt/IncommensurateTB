% Compute Chebyshev weights <v|T_n(H)|v> for use in density of states
% calculuation

addpath("chebyshev");

p = 1000;           % # polynomials in Chebyshev expansion
N = 4;              % # shifts per dimension
E_range = 13;       % energy range
r_cut = 100;        % Radial truncation
theta = 6*pi/180;   % Rotation angle
filename = 'r100_N4_p1000_dos.mat';

[X,Y] = meshgrid(0:(N-1),0:(N-1));
X = X/N-1/2;
Y = Y/N-1/2;
X = X(:);
Y = Y(:);

f = @(x) norm(x) < r_cut; % intralayer cut-off function

[S,O] = meshgrid(1:2,1:2); % fix to just 1 if you want LDOS
S = S(:);
O = O(:);
cheb_wgts = zeros(p,size(X(:),1)*4);

tstart = tic;

loc_s = graphene_init(theta,f,f,r_cut); % loc_s stores geometry, hopping function, system information.
fprintf('begin loop\n')
parfor k = 1:size(X(:),1)*4
    [i,j] = ind2sub([size(X(:),1),4],k);
    disp(['Shift = ',num2str(i),', Sheet/orbital = ',num2str(j)])

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

    disp('Generating Chebyshev weights...')
    tic;
    cheb_wgts(:,k) = Cheb_LDoS_Weights(H, E_range, v, p-1);
    disp(['Time=',num2str(toc)])
end
cheb_wgts = reshape(cheb_wgts,[p,size(X(:),1),4]);

disp(['Total time=',num2str(toc(tstart))])

save(['cheb_wgts_data/',filename], 'N', 'E_range', 'r_cut', 'theta', 'cheb_wgts');