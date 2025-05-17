% Generate AB-stacked bilayer graphene over (2*L+1)^2 lattice,
% delta = gating, gamma = interlayer coupling strength.
% fixed a connection bug on the Tx,Ty terms, Lambda --> Lambda'

function H = Generate_BLG(L, delta, gamma)

N = 2*L+1;
sigma1 = [0 1; 1 0];
Lambda = [0 1; 0 0];
sigma3 = [1 0; 0 -1];
Tx = Generate_Tx(N);
Ty = Generate_Ty(N);
I = SparseIdentity(N^2);

H = kron( delta*kron(sigma3,eye(2))+gamma*kron(Lambda,Lambda')+gamma*kron(Lambda',Lambda)...
    + kron(eye(2),sigma1) , I);
Hp = kron(kron(eye(2),Lambda'),Tx) + kron(kron(eye(2),Lambda'),Ty);
H = H + Hp + Hp';

end