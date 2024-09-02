
% for input array E and polynomial order P, calculate all T_n(E) for n <= P

function T_E = Cheb_Eval(E, P)   % E can be an array, but P is a scalar

N = size(E,2);

T_E = zeros(P+1,N);

T_E(1,:) = ones(1,N);
T_E(2,:) = E;

for j = 2:P
    T_E(j+1,:) = 2*E.*T_E(j,:) - T_E(j-1,:);
end

end