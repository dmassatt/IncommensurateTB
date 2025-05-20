% This function computes the local density of states using the Jackson kernel
% polynomial method (KPM)
%
% Inputs:
%
% p:    # polynomials in Chebyshev expansion (degree is p-1)
% oms:  Frequencies at which to evaluate the LDOS
% cheb_wgts: Chebyshev weights <v | T_n(H) | v> for all n < p. To compute
% multiple LDOS, this should be a matrix with the rows indexing n, and the
% columns indexing the different local densities of states.
%
% Note: we assume the spectrum of the Hamiltonian is on [-1,1], so the input
% frequencies have been shifted and scaled accordingly.
%
% Output:
%
% ldos: Local density of states; the rows index the frequencies oms, and the
% columns index the different local densities of states
function ldos = kpm_ldos(p, oms, E_range, cheb_wgts)

  omssc = oms/E_range; % Scale the frequencies to [-1,1]
  jackson_coeff = Cheb_JacksonCoeff(p-1);
  measure_weight = 1./sqrt(1 - omssc.^2);
  cheb_energy = Cheb_Eval(omssc, p-1);
  d = [.5 ones(1,p-1)];
  cheb_energy = diag(d)*cheb_energy;
  ldos = (((jackson_coeff.*cheb_wgts.') * cheb_energy) .*measure_weight)';
  ldos = 2*ldos / (E_range*pi);

end