% Generate AB-stacked bilayer graphene over (2*L+1)^2 lattice,
% delta = gating, gamma = interlayer coupling strength.
% Includes second nearest neighbor intralayer term
% and includes second nearest neighbor interlayer term

function H = Generate_BLG_2NN(L,delta,gamma0,gamma1,gamma3,gamma4)

% gamma0 = 3.1, gamma1 = .4, gamma3 = .3, gamma4 = .04;

% gamma0 = NN intralayer
% gamma1 = NN A-B site (on-top)
% gamma3 = NN B-A site (not on-top)
% gamma4 = 2nd NN A-A site (on-top next over)
N = 2*L+1;
sigma1 = [0 1; 1 0];
Lambda = [0 1; 0 0];
%G = [0 0; 0 1];
sigma3 = [1 0; 0 -1];
Tx = Generate_Tx(N);
Ty = Generate_Ty(N);
I = SparseIdentity(N^2);

Hgate = kron(delta*kron(sigma3,eye(2)),I);

H0 = gamma0*kron(kron(eye(2),Lambda),Tx+Ty);
H0 = H0+H0'+ gamma0*kron(kron(eye(2),sigma1) , I);

H1 = gamma1*kron(kron(Lambda,Lambda'),I);
H1 = H1 + H1';

H3 = gamma3*kron(kron(Lambda,Lambda),Tx+Ty+I);
H3 = H3 + H3';

H4 = gamma4*kron(kron(Lambda,eye(2)),Tx+Ty+I);
H4 = H4+H4';

H = Hgate+H0+H1+H3+H4;

end