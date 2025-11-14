
% inputs are the following:
% H = Hamiltonian matrix, already rescaled to have spectra in (-1,1).
% v, the "site" of interest, though it need not be local.
% P, the Chebyshev Polynomial order. Jackson coefficients will be used.
% E, an array of energies to be sampled. All E in (-1,1).

% Output will be an array LDoS, where LDoS is a function of energy
% mathematically.

function [ldos,cheb_weights] = Cheb_LDoS(H, E_range, v, P, E)
cheb_weights = Cheb_LDoS_Weights(H, E_range, v, P);
jackson_coeff = Cheb_JacksonCoeff(P+1);
measure_weight = 1./sqrt(1 - E.^2);
cheb_energy = Cheb_Eval(E, P);

d = [.5 ones(1,P)];
cheb_energy = diag(d)*cheb_energy;

ldos =  (  (jackson_coeff.*cheb_weights) * cheb_energy  ) .*measure_weight;
end