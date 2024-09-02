
% for a given matrix H, vector v, and polynomial order P, this function
% computes < v | T_n(H) | v > for all n <= P and returns it in an array of
% size P+1, coefficients n = 0, 1, . . . P.

function Cheb_Coeff = Cheb_LDoS_Weights(H, v, P)

v_1 = v;
v_2 = H*v;

Cheb_Coeff = zeros(1,P+1);
Cheb_Coeff(1) = v'*v_1;         % compute first two Chebyshev LDoS weights
Cheb_Coeff(2) = v'*v_2;

for j = 2:P
    v_new = 2*H*v_2 - v_1;
    Cheb_Coeff(j+1) = v'*v_new;
    v_1 = v_2;
    v_2 = v_new;
end

end