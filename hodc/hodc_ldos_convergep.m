% This function computes the local density of states using the high-order delta
% Chebyshev method.
%
% Given a broadening parameter eta, it performs a self-convergence test to
% determine the minimal number of Chebyshev polynomials needed to achieve a
% result converged to within a specified tolerance.
%
% Inputs:
%
% m:    Order of method with respect to broadening parameter eta
% eta:  Broadening parameter
% oms:  Frequencies at which to evaluate the LDOS
% E_range: Spectral radius of Hamiltonian
% cheb_wgts: Chebyshev weights <v | T_n(H) | v> for all n < p. To compute
% multiple LDOS, this should be a matrix with the rows indexing n, and the
% columns indexing the different local densities of states.
% tol:  Tolerance for self-convergence
% pstep: Step size for increasing the number of Chebyshev polynomials within
% self-convergence
%
% Note: we assume the spectrum of the Hamiltonian is on [-1,1], so the input
% frequencies have been shifted and scaled accordingly.
%
% Output:
%
% ldos: Local density of states; the rows index the frequencies oms, and the
% columns index the different local densities of states
% p:    Number of Chebyshev polynomials needed to achieve self-convergence; if
% self-convergence is not achieved, p will be set to the largest value tried.
function [ldos, p] = hodc_ldos_convergep(m, eta, oms, E_range, cheb_wgts, tol, pstep)

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

  pmax = size(cheb_wgts, 1);
  coefs = reshape(lorentz_coeffs(pmax,nus.',eta,-E_range,E_range),pmax,nom,m);
  for l=1:m
      coefs(:,:,l) = delta_wgt(l)*coefs(:,:,l);
  end

  % Compute the local densities of states
  nldos = size(cheb_wgts, 2);
  ldosprev = zeros(nom, nldos);

  fprintf('pstep=%d\n', pstep);
  for p=1:pstep:pmax
  ldos = zeros(nom, nldos);
    for l=1:m
      ldos = ldos + imag(coefs(1:p, :, l).' * cheb_wgts(1:p));
    end
    selferr = max(max(abs(ldos-ldosprev)));
    fprintf('p=%d, self error=%e\n', p, selferr);
    if selferr < tol
      break;
    end
    ldosprev = ldos;
  end

end