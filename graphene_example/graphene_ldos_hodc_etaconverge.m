% Compute graphene local density of states using high-order delta-Chebyshev method, for
% decreasing values of eta, and plot self-convergence error.
% Local Chebyshev weights <v|T_n(H)|v> computed using get_cheb_wgts_ldos script.
% We load a vector cheb_wgts of weights produced by that script.

% Input parameters
p = 16000;      % Chebyshev degree
m = 4;          % Order of method with respect to broadening parameter eta
dE = 0.002;    % Energy grid spacing

rate = 10^(1/m);
etas = 0.2./rate.^(0:10);   % Broadening parameters

Ls = [200 400 800 1600 3200];   % System size parameters

E = 2.0; % Energy at which to measure error

pdata = 16000;  % Chebyshev degree computed in data file

addpath('../hodc','../hodc/kernels');

[ldos_val,err] = deal(zeros(length(etas),length(Ls)));
for j=1:length(Ls)
    L = Ls(j);
    filename = ['graphene_L',num2str(L),'_p',num2str(pdata),'_ldos.mat'];
    load(['cheb_wgts_data/',filename]); % Load parameters and Chebyshev weights from file

    % Compute local densities of states
    for i=1:length(etas)
        eta = etas(i);
        ldos_val(i,j) = hodc_ldos(m, eta, p, E, E_range, cheb_wgts(1:p));

        % Compute error
        err(i,j) = abs(ldos_val(i,j) - graphene_analytic(E,1));
    end
end

% err = abs(ldos_val(2:end,:) - ldos_val(1:end-1,:));
figure(1);
set(gcf,'position',[100,100,1000,800])
loglog(etas,err,'.-','linewidth',1.5,'markersize',20); hold on
loglog(etas(2:5),1e-3*etas(2:5).^m,'--k'); hold off
xlabel('$\eta$','interpreter','latex')
ylabel('Self-convergence error')
set(gca,'fontsize',15)

% Legend with text 'r = {Ls(1),...,Ls(end)}'
legend_str = cell(1,length(Ls));
for j=1:length(Ls)
    legend_str{j} = ['L = ',num2str(Ls(j))];
end
legend(legend_str,'location','northwest')
set(gca,'fontsize',15)
