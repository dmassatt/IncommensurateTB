% Compute TBG local density of states using KPM with Jackson smoothing.
% Measure convergence with respect to p.
% Local Chebyshev weights <v|T_n(H)|v> computed using get_cheb_wgts_ldos script.
% We load a vector cheb_wgts of weights produced by that script.
%
% We generate a plot of the LDOS for different values of p, and compute
% the value of the LDOS at a specific point for different values of p.

% Input parameters
pdata = 16000;  % Chebyshev degree computed in data file

rs = [100 200 400 800 1600 3200];
ps = [500, 1000, 2000, 4000, 8000, 16000];  % Polynomial degrees
E = -0.43; % Pick specific energy to measure

addpath('../kpm')

ldos_val = zeros(length(ps),length(rs));
for j=1:length(rs)
    r = rs(j);
    filename = ['r',num2str(r),'_p',num2str(pdata),'_ldos.mat'];
    load(['../cheb_wgts_data/',filename]); % Load parameters and Chebyshev weights from file

    % Compute local densities of states
    for i=1:length(ps)
        p = ps(i);
        ldos_val(i,j) = kpm_ldos(p, E/E_range, cheb_wgts(1:p));
    end
end

% Compute self-convergence error
err = abs(ldos_val(2:end,:) - ldos_val(1:end-1,:));
figure(5);
set(gcf,'position',[100,100,1000,800])
loglog(ps(1:end-1),err,'.-','linewidth',1.5,'markersize',20); hold on
loglog(ps(2:5),1e2*ps(2:5).^-2,'--k'); hold off
xlabel('$p$','interpreter','latex')
ylabel('Self-convergence error')
set(gca,'fontsize',15)
ylim([1e-8, 1e-1])
xlim([5e2 1.5e4])

% Legend with text 'r = {rs(1),...,rs(end)}'
legend_str = cell(1,length(rs));
for j=1:length(rs)
    legend_str{j} = ['r = ',num2str(rs(j))];
end
legend(legend_str,'location','southwest')
set(gca,'fontsize',15)