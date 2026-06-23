% plot_metric_by_case.m
% Generate RMSE/MAE-vs-noise-condition figures using the same styles as generate_all_RMSE_figures.m

clear; close all; clc;

% Methods in the order of data rows (8 methods, no CRLB)
methods = {'KF', 'PF', 'Huber', 'RKF-Gaussian', 'RKF-ST', ...
           'RKF-GSTM', 'RKF-TAU', 'RKF-AGST'};

% Line style strings from generate_all_RMSE_figures (first 8, without CRLB)
line_styles = {'--^b', '--sk', '--dg', '--+m', '--xc', '-.b', '-.g', '-or'};

% Parse each style string into LineSpec components
styles = cell(1, length(line_styles));
for i = 1:length(line_styles)
    s = line_styles{i};
    % Extract line style
    if contains(s, '--')
        ls = '--';
    elseif contains(s, '-.')
        ls = '-.';
    elseif contains(s, ':')
        ls = ':';
    else
        ls = '-';
    end
    % Extract marker (second last character)
    marker = s(end-1);
    % Extract color (last character)
    color = s(end);
    styles{i} = struct('LineStyle', ls, 'Marker', marker, 'Color', color);
end

cases = {'Mild', 'Moderate', 'Severe', 'Extreme'};

% ---------- AGST data (8×4) ----------
AGST_RMSE = [
    2.1010 3.1179 4.7343 10.2260;   % KF
    2.1260 3.1556 4.8371 10.6000;   % PF
    2.2203 3.4694 5.7206 11.6880;   % Huber
    2.1009 3.1159 4.7131 10.1990;   % RKF-Gaussian
    2.1040 3.1417 4.7150 7.4064;    % RKF-ST
    2.1009 3.0898 4.5115 8.2880;    % RKF-GSTM
    2.5877 4.0657 6.8621 12.7650;   % RKF-TAU
    2.1519 3.1556 4.3826 6.3288];   % RKF-AGST

AGST_MAE = [
    1.8588 2.7077 4.0070 8.8937;
    1.8738 2.7384 4.1255 9.2330;
    1.9581 2.9656 4.5808 8.1953;
    1.8587 2.7068 3.9970 8.8814;
    1.8612 2.7247 3.9936 6.1957;
    1.8587 2.6927 3.9090 7.2588;
    2.2732 3.4305 5.3252 8.7267;
    1.8951 2.7596 3.7814 5.4180];

% ---------- GMM data (8×4) ----------
GMM_RMSE = [
    1.310 2.363 3.820 4.419;
    1.334 2.448 3.880 4.531;
    1.352 2.539 4.399 5.238;
    1.310 2.360 3.813 4.407;
    1.314 2.539 4.371 4.597;
    1.310 2.239 3.546 3.953;
    1.515 2.988 5.156 6.335;
    1.316 2.305 3.116 3.442];

GMM_MAE = [
    1.147 1.906 3.123 3.626;
    1.164 1.986 3.197 3.755;
    1.177 1.898 3.150 3.824;
    1.147 1.906 3.119 3.619;
    1.149 1.908 3.167 3.554;
    1.147 1.825 2.910 3.293;
    1.313 2.083 3.311 4.268;
    1.155 1.779 2.359 2.851];

% ---------- Laplace Outliers data (8×4) ----------
Laplace_RMSE = [
    2.377 5.154 13.376 24.677;
    2.575 5.452 13.760 24.870;
    2.490 6.007 21.177 53.300;
    2.359 5.082 13.186 24.338;
    2.117 2.855 4.567 7.144;
    1.650 2.802 7.198 14.142;
    2.934 7.041 24.363 62.517;
    1.599 1.810 3.035 4.364];

Laplace_MAE = [
    1.734 3.908 11.012 20.506;
    1.962 4.262 11.449 20.752;
    1.543 3.148 11.310 27.513;
    1.729 3.876 10.894 20.277;
    1.489 1.850 2.770 4.336;
    1.405 2.377 6.204 12.238;
    1.716 3.178 10.576 28.421;
    1.293 1.410 2.080 2.877];

% ---------- Plot function ----------
function plot_one(Y, methods, styles, cases, ylab, ttl, fname)
    figure('Color','w','Position',[100,100,780,480]);
    hold on; grid on;
    x = 1:numel(cases);
    for i = 1:numel(methods)
        s = styles{i};
        plot(x, Y(i,:), 'LineStyle', s.LineStyle, 'Marker', s.Marker, ...
             'Color', s.Color, 'LineWidth', 1.4, 'MarkerSize', 6);
    end
    set(gca, 'XTick', x, 'XTickLabel', cases);
    xlabel('Noise condition');
    ylabel(ylab);
    title(ttl);
    legend(methods, 'Location', 'northwest');
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 11);
    exportgraphics(gcf, [fname '.png'], 'Resolution', 300);
    print(gcf, [fname '.eps'], '-depsc');
    close(gcf);
end

% Generate all figures
plot_one(AGST_RMSE, methods, styles, cases, 'RMSE (m)', ...
         'RMSE under AGST noises', 'AGST_RMSE_by_case');
plot_one(AGST_MAE, methods, styles, cases, 'MAE (m)', ...
         'MAE under AGST noises', 'AGST_MAE_by_case');
plot_one(GMM_RMSE, methods, styles, cases, 'RMSE (m)', ...
         'RMSE under GMM noises', 'GMM_RMSE_by_case');
plot_one(GMM_MAE, methods, styles, cases, 'MAE (m)', ...
         'MAE under GMM noises', 'GMM_MAE_by_case');
plot_one(Laplace_RMSE, methods, styles, cases, 'RMSE (m)', ...
         'RMSE under Gaussian + Laplace outliers', 'Laplace_RMSE_by_case');
plot_one(Laplace_MAE, methods, styles, cases, 'MAE (m)', ...
         'MAE under Gaussian + Laplace outliers', 'Laplace_MAE_by_case');

fprintf('All figures generated.\n');