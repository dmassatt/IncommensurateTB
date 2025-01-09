% Load data files
r_values = [100, 200, 400, 800, 1600];
E_data = cell(1, length(r_values));
for idx = 1:length(r_values)
    filename = sprintf('E0_ldos_values_kpm_r%d.dat', r_values(idx));
    E_data{idx} = load(filename);
end

% Compute differences
diffs = cell(1, length(r_values));
for idx = 1:length(r_values)
    E = E_data{idx};
    diffs{idx} = abs(E(1:end-1, 2) - E(2:end, 2));
end

% Create log-log plot
figure('Position', [100, 100, 800, 600]);
colors = lines(length(r_values)); % Generate distinct colors for each line
for idx = 1:length(r_values)
    E = E_data{idx};
    loglog(E(1:end-1, 1), diffs{idx}, 'o-', 'DisplayName', sprintf('r=%d, p=8000', r_values(idx)), 'Color', colors(idx, :),'LineWidth', 2);
    hold on;
end

% Add slopes
E800 = E_data{4}; % r=800
incept = log10(diffs{4}(2)) + 2 * log10(E800(2, 1));
N_lst = E800(1:6, 1);
y = 10.^(-2 * log10(N_lst) + incept);
loglog(N_lst, y, 'k:', 'LineWidth', 2, 'DisplayName', 'slope=2');

incept = log10(diffs{4}(2)) + 3 * log10(E800(2, 1));
y = 10.^(-3 * log10(N_lst) + incept);
loglog(N_lst, y, 'b:', 'LineWidth', 2, 'DisplayName', 'slope=3');

% Set labels and legend
xlabel('$N_c$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('Self-Convergence Error', 'FontSize', 14);
set(gca, 'FontSize', 12);
legend('Location', 'best');

% Optional: save the plot
saveas(gcf, 'TBG_kpm_error_eta_p8000.png');
