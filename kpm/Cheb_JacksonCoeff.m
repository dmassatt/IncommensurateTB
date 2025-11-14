
% generates the jackson coefficients from 0 to P

function jackson = Cheb_JacksonCoeff(P)

p = 0:P-1;

jackson = ( (P - p + 1).* cos(pi*p/(P+1)) + sin(pi*p/(P+1))*cot(pi/(P+1)) ) / (P+1);


end