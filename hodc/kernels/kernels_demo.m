%% Experiment 1 - multiplication operator with continuous spectrum
X=-2:0.02:2;                               	% evaluation pts
f=chebfun(@(x) sqrt(3/2)*x);               	% measure wrt f(x)
m=chebfun(@(x) x);                          % "multiplication by" operator

% Compute spectral measures
eta = 0.001;                                % smoothing parameter
ord = 4;                                    % smoothing kernel order
[pol,res]=rational_kernel(ord,'equi');      % poles/residues of kernel
mu = zeros(size(X));
for j = 1:length(X)
    for k = 1:ord
        z = X(j) - eta*pol(k);                  % shift
        u = f / (m-z);                          % shifted linear solve
        mu(j) = mu(j) - res(k)*imag(f'*u)/pi;   % update linear combo
    end
end

% Plot
figure(1)
semilogy(X,mu,'LineWidth',2)
hold on
xlim([X(1) X(end)])
ax = gca; ax.FontSize = 14;

%% Experiment 2 - differential operator with discrete spectrum
X=0:0.1:15;                                 % evaluation pts
f=chebfun(@(x) x + x.^2);               	% measure wrt f(x)
H=chebop(@(u) -diff(u,2));                  % differential operator
H.bc = 0;                                   % boundary conditions

% Compute spectral measures
eta = 0.01;                                 % smoothing parameter
ord = 2;                                    % smoothing kernel order
[pol,res]=rational_kernel(ord,'equi');      % poles/residues of kernel
I = chebop(@(u) u);                        	% identity, for shifts
mu = zeros(size(X));
for j = 1:length(X)
    for k = 1:ord
        z = X(j) - eta*pol(k);                  % shift
        u = (H - z*I) \ f;                      % shifted linear solve
        mu(j) = mu(j) - res(k)*imag(f'*u)/pi;   % update linear combo
    end
end

% Plot
figure(2)
semilogy(X,mu,'LineWidth',2)
hold on
xlim([X(1) X(end)])
ax = gca; ax.FontSize = 14;