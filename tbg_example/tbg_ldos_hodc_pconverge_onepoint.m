% Compute TBG local density of states using KPM with Jackson smoothing, for
% decreasing values of eta, and plot self-convergence error at a single energy point.
% Local Chebyshev weights <v|T_n(H)|v> computed using get_cheb_wgts_ldos script.
% We load a vector cheb_wgts of weights produced by that script.

% Input parameters
m = 6;          % Order of method with respect to broadening parameter eta
pdata = 16000;  % Chebyshev degree computed in data file

rate = 10^(1/m);
rs = [100 200 400 800 1600 3200];
etas = 0.02./rate.^(0:9);   % Broadening parameters
E = -0.43; % Pick specific energy to measure

tol = 1e-8;
pstep = 50;

addpath('../hodc','../hodc/kernels');

[ldos_val,p] = deal(zeros(length(etas),length(rs)));
for j=1:length(rs)
    r = rs(j);
    filename = ['r',num2str(r),'_p',num2str(pdata),'_ldos.mat'];
    load(['../cheb_wgts_data/',filename]); % Load parameters and Chebyshev weights from file

    % Compute local densities of states
    for i=1:length(etas)
        eta = etas(i);
        fprintf('Computing LDOS for r = %d, eta = %f\n', r, eta);
        [ldos_val(i,j),p(i,j)] = hodc_ldos_convergep(m, eta, E/E_range, cheb_wgts(1:pdata), tol, pstep);
    end
end

% Compute self-convergence error
err = abs(ldos_val(2:end,:) - ldos_val(1:end-1,:));
figure(3);
set(gcf,'position',[100,100,1000,800])
for j=1:length(rs)
    loglog(p(1:end-1,j),err(:,j),'.-','linewidth',1.5,'markersize',20); hold on
end
loglog(p(3:6,end),1e15*p(3:6,end).^(-m),'--k');
loglog(p(3:6,end),1e2*p(3:6,end).^(-2),'--k');
hold off
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
