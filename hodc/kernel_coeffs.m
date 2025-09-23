% Chebyshev coefficients for the kernel function
%
% K(x, omega) = -1/pi * sum_{l=1}^m Im( w_l / (omega - x + eta * x_l + i * eta)^-1 )
%
% on [a,b]
function [coefs] = kernel_coeffs(p, om, eta, delta_polx, delta_wgt, a, b)

  pp = 4*p; % Oversampling factor

  % Get Chebyshev nodes on [a, b]
  xc = cos(pi*(2*(0:pp-1)' + 1)/(2*pp)); % Chebyshev nodes on [-1, 1]
  xc = (a + b)/2 + (b - a)/2*xc; % Chebyshev nodes on [a, b]

  % Evaluate the Lorentzian function at the Chebyshev nodes
  w = reshape(delta_wgt, 1, 1, []);
  x = reshape(delta_polx, 1, 1, []);
  f = -1 / pi * imag(sum(w ./ (om - xc + (eta * x + 1i * eta)), 3));

  % Compute the Chebyshev coefficients
  coefs = dct(f,'Type',2);
  coefs(1,:) = coefs(1,:)/sqrt(2);

  coefs = sqrt(2/pp)*coefs(1:p,:);

end