% Chebyshev coefficients for the Lorentzian function
%
% f(x) = -1/pi * (nu - x + i*eta)^-1
%
% on [a,b]
function [coefs] = lorentz_coeffs_vec(p, nu, eta, a, b)

  pp = 4*p; % Oversampling factor

  % Get Chebyshev nodes on [a, b]
  xc = cos(pi*(2*(0:pp-1)' + 1)/(2*pp)); % Chebyshev nodes on [-1, 1]
  xc = (a + b)/2 + (b - a)/2*xc; % Chebyshev nodes on [a, b]

  % Evaluate the Lorentzian function at the Chebyshev nodes
  f = -1/pi./(nu - xc + 1i*eta);

  % Compute the Chebyshev coefficients
  coefs = sqrt(2/pp)*dct(f,'Type',2);
  coefs(1,:) = coefs(1,:)/sqrt(2);

  coefs = coefs(1:p,:);

end