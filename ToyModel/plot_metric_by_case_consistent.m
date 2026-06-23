% plot_metric_by_case_consistent.m
% Generate RMSE/MAE-vs-noise-condition figures with a unified thesis style.
% Main change: use one global grid switch for all figures.
% Recommended for the current thesis: USE_GRID = false, so these figures
% match the cleaner no-grid style of the RMSE evolution figures.

clear; close all; clc;

% ===================== Global visual style =====================
USE_GRID = false;                 % false: no grid; true: subtle grid for all figures
FIG_POS  = [100, 100, 760, 460];  % consistent aspect ratio
FONT_NAME = 'Times New Roman';
FONT_SIZE = 11;
AX_LINE_WIDTH = 0.8;
BASE_LINE_WIDTH = 1.25;
OURS_LINE_WIDTH = 1.65;
MARKER_SIZE = 6;

% Methods in the order of data rows (8 methods, no CRLB)
methods = {'KF', 'PF', 'Huber', 'RKF-Gaussian', 'RKF-ST', ...
           'RKF-GSTM', 'RKF-TAU', 'RKF-AGST'};

% Styles follow the old generate_all_RMSE_figures.m convention.
% Note: '-.b' and '-.g' mean dash-dot lines without markers.
style_specs = {'--^b', '--sk', '--dg', '--+m', '--xc', '-.b', '-.g', '-or'};
styles = cell(1, numel(style_specs));
for i = 1:numel(style_specs)
    styles{i} = parse_line_spec(style_specs{i});
end

cases = {'Mild', 'Moderate', 'Severe', 'Extreme'};

% ===================== AGST data (8 x 4) =====================
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

% ===================== GMM data (8 x 4) =====================
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

% ===================== Gaussian + Laplace outlier data (8 x 4) =====================
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

% ===================== Generate all figures =====================
plot_one(AGST_RMSE, methods, styles, cases, 'RMSE (m)', ...
         'RMSE under AGST noises', 'AGST_RMSE_by_case', ...
         USE_GRID, FIG_POS, FONT_NAME, FONT_SIZE, AX_LINE_WIDTH, ...
         BASE_LINE_WIDTH, OURS_LINE_WIDTH, MARKER_SIZE);

plot_one(AGST_MAE, methods, styles, cases, 'MAE (m)', ...
         'MAE under AGST noises', 'AGST_MAE_by_case', ...
         USE_GRID, FIG_POS, FONT_NAME, FONT_SIZE, AX_LINE_WIDTH, ...
         BASE_LINE_WIDTH, OURS_LINE_WIDTH, MARKER_SIZE);

plot_one(GMM_RMSE, methods, styles, cases, 'RMSE (m)', ...
         'RMSE under GMM noises', 'GMM_RMSE_by_case', ...
         USE_GRID, FIG_POS, FONT_NAME, FONT_SIZE, AX_LINE_WIDTH, ...
         BASE_LINE_WIDTH, OURS_LINE_WIDTH, MARKER_SIZE);

plot_one(GMM_MAE, methods, styles, cases, 'MAE (m)', ...
         'MAE under GMM noises', 'GMM_MAE_by_case', ...
         USE_GRID, FIG_POS, FONT_NAME, FONT_SIZE, AX_LINE_WIDTH, ...
         BASE_LINE_WIDTH, OURS_LINE_WIDTH, MARKER_SIZE);

plot_one(Laplace_RMSE, methods, styles, cases, 'RMSE (m)', ...
         'RMSE under Gaussian + Laplace outliers', 'Laplace_RMSE_by_case', ...
         USE_GRID, FIG_POS, FONT_NAME, FONT_SIZE, AX_LINE_WIDTH, ...
         BASE_LINE_WIDTH, OURS_LINE_WIDTH, MARKER_SIZE);

plot_one(Laplace_MAE, methods, styles, cases, 'MAE (m)', ...
         'MAE under Gaussian + Laplace outliers', 'Laplace_MAE_by_case', ...
         USE_GRID, FIG_POS, FONT_NAME, FONT_SIZE, AX_LINE_WIDTH, ...
         BASE_LINE_WIDTH, OURS_LINE_WIDTH, MARKER_SIZE);

fprintf('All figures generated with USE_GRID = %d.\n', USE_GRID);

% ===================== Local functions =====================
function plot_one(Y, methods, styles, cases, ylab, ttl, fname, ...
                  useGrid, figPos, fontName, fontSize, axLineWidth, ...
                  baseLineWidth, oursLineWidth, markerSize)
    fig = figure('Color','w','Position',figPos, 'Renderer','painters');
    ax = axes(fig);
    hold(ax, 'on');

    x = 1:numel(cases);
    for i = 1:numel(methods)
        s = styles{i};
        lw = baseLineWidth;
        if strcmp(methods{i}, 'RKF-AGST')
            lw = oursLineWidth;
        end
        plot(ax, x, Y(i,:), ...
             'LineStyle', s.LineStyle, ...
             'Marker', s.Marker, ...
             'Color', s.Color, ...
             'LineWidth', lw, ...
             'MarkerSize', markerSize, ...
             'MarkerFaceColor', 'none');
    end

    set(ax, 'XTick', x, 'XTickLabel', cases);
    xlim(ax, [0.85, numel(cases) + 0.15]);
    xlabel(ax, 'Noise condition');
    ylabel(ax, ylab);
    title(ax, ttl, 'FontWeight', 'normal');

    format_axes(ax, useGrid, fontName, fontSize, axLineWidth);

    lgd = legend(ax, methods, 'Location', 'northwest');
    set(lgd, 'Box', 'on', 'FontName', fontName, 'FontSize', fontSize - 1);

    % Add a small y-axis margin without changing the overall scale too much.
    yl = ylim(ax);
    dy = yl(2) - yl(1);
    if dy > 0
        ylim(ax, [max(0, yl(1) - 0.03 * dy), yl(2) + 0.08 * dy]);
    end

    % Export both raster and vector versions.
    exportgraphics(fig, [fname '.png'], 'Resolution', 300);
    print(fig, [fname '.eps'], '-depsc', '-painters', '-r300');
    close(fig);
end

function format_axes(ax, useGrid, fontName, fontSize, axLineWidth)
    set(ax, ...
        'FontName', fontName, ...
        'FontSize', fontSize, ...
        'LineWidth', axLineWidth, ...
        'Box', 'on', ...
        'TickDir', 'in', ...
        'XMinorTick', 'off', ...
        'YMinorTick', 'off', ...
        'XMinorGrid', 'off', ...
        'YMinorGrid', 'off', ...
        'Layer', 'top');

    if useGrid
        grid(ax, 'on');
        set(ax, 'GridLineStyle', ':', 'GridAlpha', 0.20);
    else
        grid(ax, 'off');
    end
end

function st = parse_line_spec(spec)
    % Robust parser for simple MATLAB LineSpec strings such as '--^b', '-.b', '-or'.
    colors = 'bgrcmykw';
    markers = {'+', 'o', '*', '.', 'x', 's', 'd', '^', 'v', '>', '<', 'p', 'h'};

    if contains(spec, '--')
        lineStyle = '--';
        rest = erase(spec, '--');
    elseif contains(spec, '-.')
        lineStyle = '-.';
        rest = erase(spec, '-.');
    elseif contains(spec, ':')
        lineStyle = ':';
        rest = erase(spec, ':');
    elseif contains(spec, '-')
        lineStyle = '-';
        rest = erase(spec, '-');
    else
        lineStyle = '-';
        rest = spec;
    end

    color = 'k';
    marker = 'none';
    for j = 1:length(rest)
        ch = rest(j);
        if contains(colors, ch)
            color = ch;
        elseif any(strcmp(ch, markers))
            marker = ch;
        end
    end

    st = struct('LineStyle', lineStyle, 'Marker', marker, 'Color', color);
end
