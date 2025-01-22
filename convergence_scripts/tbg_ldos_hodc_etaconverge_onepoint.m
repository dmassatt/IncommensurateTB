% Compute TBG local density of states using KPM with Jackson smoothing, for
% decreasing values of eta, and plot self-convergence error at a single energy point.
% Local Chebyshev weights <v|T_n(H)|v> computed using get_cheb_wgts_ldos script.
% We load a vector cheb_wgts of weights produced by that script.

% Input parameters
p = 16000;       % Chebyshev degree
m = 6;          % Order of method with respect to broadening parameter eta
pdata = 16000;  % Chebyshev degree computed in data file

rate = 10^(1/m);
rs = [100 200 400 800 1600 3200];
etas = 0.02./rate.^(0:10);   % Broadening parameters
E = -0.43; % Pick specific energy to measure

addpath('../hodc','../hodc/kernels');

ldos_val = zeros(length(etas),length(rs));
for j=1:length(rs)
    r = rs(j);
    filename = ['r',num2str(r),'_p',num2str(pdata),'_ldos.mat'];
    load(['../cheb_wgts_data/',filename]); % Load parameters and Chebyshev weights from file

    % Compute local densities of states
    for i=1:length(etas)
        eta = etas(i);
        ldos_val(i,j) = hodc_ldos(m, eta, p, E/E_range, cheb_wgts(1:p));
    end
end

% Compute self-convergence error
err = abs(ldos_val(2:end,:) - ldos_val(1:end-1,:));
figure(2);
loglog(etas(1:end-1),err,'.-','linewidth',1.5,'markersize',20); hold on
loglog(etas(2:5),1e10*etas(2:5).^m,'--k'); hold off
xlabel('$\eta$','interpreter','latex')
ylabel('Self-convergence error')
set(gca,'fontsize',15)

% Legend with text 'r = {rs(1),...,rs(end)}'
legend_str = cell(1,length(rs));
for j=1:length(rs)
    legend_str{j} = ['r = ',num2str(rs(j))];
end
legend(legend_str,'location','northwest')
set(gca,'fontsize',15)
