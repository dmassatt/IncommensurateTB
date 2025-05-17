% Generate AB-stacked bilayer graphene over (2*L+1)^2 lattice,
% delta = gating, gamma = interlayer coupling strength.
% fixed a connection bug on the Tx,Ty terms, Lambda --> Lambda'

function H = Generate_MonG(L, t)

N = 2*L+1;
sigma1 = [0 1; 1 0];
Lambda = [0 1; 0 0];
Tx = Generate_Tx(N);
Ty = Generate_Ty(N);
I = SparseIdentity(N^2);

H = kron(  t*sigma1 , I);
Hp = kron(Lambda',Tx+Ty);
H = H + Hp + Hp';

end