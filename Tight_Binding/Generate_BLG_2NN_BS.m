function E = Generate_BLG_2NN_BS(delta,gamma0,gamma1,gamma3,gamma4, Kx, Ky)

% gamma0 = 3.1, gamma1 = .4, gamma3 = .3, gamma4 = .04;

% gamma0 = NN intralayer
% gamma1 = NN A-B site (on-top)
% gamma3 = NN B-A site (not on-top)
% gamma4 = 2nd NN A-A site (on-top next over)

N_K = size(Kx,2);
E = zeros(4,N_K);

sigma1 = [0 1; 1 0];
Lambda = [0 1; 0 0];
G = eye(2);
sigma3 = [1 0; 0 -1];
for j = 1:N_K
    Tx = exp(1i*Kx(j));
    Ty = exp(1i*Ky(j));
    
    Hgate = delta*kron(sigma3,eye(2));
    
    H0 = gamma0*kron(eye(2),Lambda)*(Tx+Ty);
    H0 = H0+H0'+ gamma0*kron(eye(2),sigma1);
    
    H1 = gamma1*kron(Lambda,Lambda');
    H1 = H1 + H1';
    
    H3 = gamma3*kron(Lambda,Lambda)*(Tx+Ty+1);
    H3 = H3 + H3';
    
    H4 = gamma4*kron(Lambda,G)*(Tx+Ty+1);
    H4 = H4+H4';
    
    H = Hgate+H0+H1+H3+H4;
    E(:,j) = eig(H);

end
end