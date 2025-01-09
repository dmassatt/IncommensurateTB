% Load data files
r_values = [100, 200, 400, 800, 1600];
E_data = cell(1, length(r_values));
for idx = 1:length(r_values)
    filename = sprintf('E0_ldos_values_hodc_r%d_p8000.dat', r_values(idx));
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
incept = log10(diffs{4}(7)) - 6 * log10(E800(7, 1));
eta_lst = E800(5:11, 1);
y = 10.^(6 * log10(eta_lst) + incept);
loglog(eta_lst, y, 'k:', 'LineWidth', 2, 'DisplayName', 'slope=6');

incept = log10(diffs{4}(7)) - 4 * log10(E800(7, 1));
y = 10.^(4 * log10(eta_lst) + incept);
loglog(eta_lst, y, 'b:', 'LineWidth', 2, 'DisplayName', 'slope=4');

incept = log10(diffs{4}(7)) - 3 * log10(E800(7, 1));
y = 10.^(3 * log10(eta_lst) + incept);
loglog(eta_lst, y, 'Color', [1, 0.5, 0], 'LineStyle', ':', 'LineWidth', 2, 'DisplayName', 'slope=3');

% Set labels and legend
xlabel('\eta', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('Self-Convergence Error', 'FontSize', 14);
set(gca, 'FontSize', 12);
legend('Location', 'best');

% Optional: save the plot
% saveas(gcf, 'TBG_hodc_error_eta_p8000.png');
