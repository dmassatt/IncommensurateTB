% This function computes the local density of states using the high-order delta
% Chebyshev method
%
% Inputs:
%
% m:    Order of method with respect to broadening parameter eta
% eta:  Broadening parameter
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
function ldos = hodc_ldos(m, eta, p, oms, E_range, cheb_wgts)

  % Get poles and weights for high-order expansion of delta function
  if (m > 6)
    warning('Poles & residues may not be accurate for m > 6.');
  end
  [delta_pol,delta_wgt]=rational_kernel(m,'equi');
  delta_polx = real(delta_pol);

  % Get weighted Lorentzian Chebyshev coefficients for all frequencies
  nom = length(oms);

  nus = zeros(nom*m,1);
  for l = 1:m
      for n = 1:nom
          nus((l-1)*nom+n) = oms(n) - eta*delta_polx(l);
      end
  end

  coefs = reshape(lorentz_coeffs(p,nus.',eta,-E_range,E_range),p,nom,m);
  for l=1:m
      coefs(:,:,l) = delta_wgt(l)*coefs(:,:,l);
  end

  % Compute the local densities of states
  tic;
  nldos = size(cheb_wgts, 2);
  ldos = zeros(nom, nldos);
  for l=1:m
    ldos = ldos + imag(coefs(:, :, l).' * cheb_wgts);
  end

end