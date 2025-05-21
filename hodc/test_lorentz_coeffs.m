% Test the lorentz_coeffs function

eta = 0.1;
a = -1;
b = 2;

%rng(1783); % Fix random seed
x = a + (b-a)*rand(100,1);
nu = 0.3;

ps = 20:20:600;
err = zeros(size(ps));
for i=1:length(ps)
  p = ps(i);
  coefs = lorentz_coeffs(p, nu, eta, a, b);

  % Evaluate Chebyshev expansion
  xsc = 2*(x - a)/(b - a) - 1; % Evaluation points on [-1,1]
  ftest = cos(acos(xsc)*(0:p-1)) * coefs;
  ftrue = -1/pi * 1./(x - nu + 1i*eta);

  % Check the error
  err(i) = max(abs(ftest - ftrue));
  disp(['Error for p = ', num2str(p), ': ', num2str(err(i))])
end

figure(1);
semilogy(ps, err, '-o')
xlabel('p')
ylabel('max error')
% addpath ~/Documents/MATLAB/export_fig/
% export_fig('errs.pdf');
