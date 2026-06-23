% plot_metric_by_case_param_labels.m
% Generate RMSE/MAE-vs-noise-parameter figures.
%
% This version follows the time-evolution figure style:
%   - grid off
%   - tick marks inward
%   - boxed axes
%   - concrete x-axis noise parameters instead of Mild/Moderate/Severe
%
% To make sure MATLAB can display the x-axis labels correctly, all tick labels
% use plain ASCII text, not Greek symbols or LaTeX commands. Multi-line x-axis
% labels are drawn manually by text(), which is more reliable than long
% XTickLabel strings when exporting EPS/PDF figures.

clear; close all; clc;

% ===================== Global visual style =====================
USE_GRID = false;                  % match time-evolution figures: no grid
FIG_POS  = [100, 100, 880, 520];   % wider figure for parameter labels
FONT_NAME = 'Times New Roman';
FONT_SIZE = 11;
XTICK_FONT_SIZE = 9;
AX_LINE_WIDTH = 0.5;               % MATLAB-like thin axis line
BASE_LINE_WIDTH = 1.25;
OURS_LINE_WIDTH = 1.65;
MARKER_SIZE = 6;

% If true, x tick labels are manually drawn as multi-line text. This is safer
% for long parameter labels and EPS export.
USE_MANUAL_XTICK_LABELS = true;

% Methods in the order of data rows: 8 methods, no CRLB.
methods = {'KF', 'PF', 'Huber', 'RKF-Gaussian', 'RKF-ST', ...
           'RKF-GSTM', 'RKF-TAU', 'RKF-AGST', 'RKF-AGSlash'};

% Styles follow the old generate_all_RMSE_figures.m convention.
% Note: '-.b' and '-.g' mean dash-dot lines without markers.
style_specs = {'--^b', '--sk', '--dg', '--+m', '--xc', '-.b', '-.g', '-or', '-sm'};
styles = cell(1, numel(style_specs));
for i = 1:numel(style_specs)
    styles{i} = parse_line_spec(style_specs{i});
end

% Numeric x positions.
xCases = 1:4;

% ===================== X-axis labels: concrete noise parameters =====================
% AGST: v_k = v_G,k + v_t,k, v_G,k ~ N(0, Rg I2), v_t,k ~ ST(0, sigma_t I2, nu_t)
AGST_xlabels = { ...
    {'R=1', 'v=10', 'σ=2'}, ...
    {'R=1', 'v=5',  'σ=5'}, ...
    {'R=1', 'v=3',  'σ=10'}, ...
    {'R=1', 'v=2',  'σ=20'}  ...
};

% GMM: v_k ~ w1 N(0, sigma1^2 I2) + w2 N(0, sigma2^2 I2)
% Plain ASCII label uses sigma2=[sigma1^2,sigma2^2].
GMM_xlabels = { ...
    {'w=[0.95,0.05]', 'σ2=[1,5]'}, ...
    {'w=[0.80,0.20]', 'σ2=[1,20]'}, ...
    {'w=[0.70,0.30]', 'σ2=[1,50]'}, ...
    {'w=[0.80,0.20]', 'σ2=[4,100]'} ...
};

% Gaussian background + sparse Laplace outliers:
% epsilon_k ~ N(0, Rbase I2), s_k ~ Bernoulli(pout), outlier scale set by sigma_out.
Laplace_xlabels = { ...
    {'R=1', 'p=0.05', 'σ=10'}, ...
    {'R=1', 'p=0.10', 'σ=20'}, ...
    {'R=1', 'p=0.20', 'σ=50'}, ...
    {'R=1', 'p=0.30', 'σ=100'} ...
};

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
plot_one(AGST_RMSE, methods, styles, xCases, AGST_xlabels, 'RMSE (m)', ...
         'RMSE under AGST noises', 'AGST_RMSE_by_case', ...
         USE_GRID, FIG_POS, FONT_NAME, FONT_SIZE, XTICK_FONT_SIZE, AX_LINE_WIDTH, ...
         BASE_LINE_WIDTH, OURS_LINE_WIDTH, MARKER_SIZE, USE_MANUAL_XTICK_LABELS);

plot_one(AGST_MAE, methods, styles, xCases, AGST_xlabels, 'MAE (m)', ...
         'MAE under AGST noises', 'AGST_MAE_by_case', ...
         USE_GRID, FIG_POS, FONT_NAME, FONT_SIZE, XTICK_FONT_SIZE, AX_LINE_WIDTH, ...
         BASE_LINE_WIDTH, OURS_LINE_WIDTH, MARKER_SIZE, USE_MANUAL_XTICK_LABELS);

plot_one(GMM_RMSE, methods, styles, xCases, GMM_xlabels, 'RMSE (m)', ...
         'RMSE under GMM noises', 'GMM_RMSE_by_case', ...
         USE_GRID, FIG_POS, FONT_NAME, FONT_SIZE, XTICK_FONT_SIZE, AX_LINE_WIDTH, ...
         BASE_LINE_WIDTH, OURS_LINE_WIDTH, MARKER_SIZE, USE_MANUAL_XTICK_LABELS);

plot_one(GMM_MAE, methods, styles, xCases, GMM_xlabels, 'MAE (m)', ...
         'MAE under GMM noises', 'GMM_MAE_by_case', ...
         USE_GRID, FIG_POS, FONT_NAME, FONT_SIZE, XTICK_FONT_SIZE, AX_LINE_WIDTH, ...
         BASE_LINE_WIDTH, OURS_LINE_WIDTH, MARKER_SIZE, USE_MANUAL_XTICK_LABELS);

plot_one(Laplace_RMSE, methods, styles, xCases, Laplace_xlabels, 'RMSE (m)', ...
         'RMSE under Gaussian + Laplace outliers', 'Laplace_RMSE_by_case', ...
         USE_GRID, FIG_POS, FONT_NAME, FONT_SIZE, XTICK_FONT_SIZE, AX_LINE_WIDTH, ...
         BASE_LINE_WIDTH, OURS_LINE_WIDTH, MARKER_SIZE, USE_MANUAL_XTICK_LABELS);

plot_one(Laplace_MAE, methods, styles, xCases, Laplace_xlabels, 'MAE (m)', ...
         'MAE under Gaussian + Laplace outliers', 'Laplace_MAE_by_case', ...
         USE_GRID, FIG_POS, FONT_NAME, FONT_SIZE, XTICK_FONT_SIZE, AX_LINE_WIDTH, ...
         BASE_LINE_WIDTH, OURS_LINE_WIDTH, MARKER_SIZE, USE_MANUAL_XTICK_LABELS);

fprintf('All figures generated with concrete x-axis noise parameters.\n');

% ===================== Local functions =====================
function plot_one(Y, methods, styles, x, xLabelLines, ylab, ttl, fname, ...
                  useGrid, figPos, fontName, fontSize, xTickFontSize, axLineWidth, ...
                  baseLineWidth, oursLineWidth, markerSize, useManualXTickLabels)
    fig = figure('Color','w','Position',figPos, 'Renderer','painters');
    ax = axes(fig);
    hold(ax, 'on');

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

    xlim(ax, [0.85, numel(x) + 0.15]);
    ylabel(ax, ylab);
    title(ax, ttl);

    format_axes(ax, useGrid, fontName, fontSize, axLineWidth);

    % Add a small y-axis margin without changing the overall scale too much.
    yl = ylim(ax);
    dy = yl(2) - yl(1);
    if dy > 0
        ylim(ax, [max(0, yl(1) - 0.03 * dy), yl(2) + 0.08 * dy]);
    end

    % Put legend after y-limit adjustment. Location follows previous figures.
    lgd = legend(ax, methods, 'Location', 'northwest');
    set(lgd, 'Box', 'on', 'FontName', fontName, 'FontSize', fontSize - 1);

    % X-axis parameter labels. Manual labels avoid problems with long labels,
    % underscores, Greek letters, and EPS export.
    set(ax, 'XTick', x);
    if useManualXTickLabels
        set(ax, 'XTickLabel', repmat({''}, 1, numel(x)));
        % Leave bottom space for 2-3 lines of parameter text.
        set(ax, 'Position', [0.11, 0.27, 0.84, 0.64]);
        add_multiline_xtick_labels(ax, x, xLabelLines, fontName, xTickFontSize);
        hxl = xlabel(ax, 'Noise parameters');
        set(hxl, 'Units', 'normalized', 'Position', [0.5, -0.23, 0]);
    else
        % Fallback: display one-line labels directly. This is less compact.
        set(ax, 'XTickLabel', flatten_label_lines(xLabelLines), 'TickLabelInterpreter', 'none');
        xtickangle(ax, 25);
        xlabel(ax, 'Noise parameters');
    end

    % Export both raster and vector versions. Use a fallback for old MATLAB.
    set(fig, 'PaperPositionMode', 'auto');
    if exist('exportgraphics', 'file') == 2
        exportgraphics(fig, [fname '.png'], 'Resolution', 300);
    else
        print(fig, [fname '.png'], '-dpng', '-r300');
    end
    print(fig, [fname '.eps'], '-depsc', '-painters', '-r300');
    close(fig);
end

function format_axes(ax, useGrid, fontName, fontSize, axLineWidth)
    % Match the existing RMSE time-evolution figures:
    %   - TickDir = 'in'
    %   - Box = 'on'
    %   - Grid = 'off' by default
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
        'TickLabelInterpreter', 'none', ...
        'Layer', 'top');

    if useGrid
        grid(ax, 'on');
        set(ax, 'GridLineStyle', ':', 'GridAlpha', 0.20);
    else
        grid(ax, 'off');
    end
end

function add_multiline_xtick_labels(ax, xVals, labelLines, fontName, fontSize)
    yl = ylim(ax);
    dy = yl(2) - yl(1);
    yText = yl(1) - 0.075 * dy;

    for ii = 1:numel(xVals)
        text(ax, xVals(ii), yText, labelLines{ii}, ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'top', ...
            'FontName', fontName, ...
            'FontSize', fontSize, ...
            'Interpreter', 'none', ...
            'Clipping', 'off');
    end
end

function flatLabels = flatten_label_lines(labelLines)
    flatLabels = cell(1, numel(labelLines));
    for ii = 1:numel(labelLines)
        thisLines = labelLines{ii};
        label = thisLines{1};
        for jj = 2:numel(thisLines)
            label = [label, ', ', thisLines{jj}]; %#ok<AGROW>
        end
        flatLabels{ii} = label;
    end
end

function st = parse_line_spec(spec)
    % Robust parser for simple MATLAB LineSpec strings such as '--^b', '-.b', '-or'.
    colors = 'bgrcmykw';
    markers = {'+', 'o', '*', '.', 'x', 's', 'd', '^', 'v', '>', '<', 'p', 'h'};

    if ~isempty(strfind(spec, '--'))
        lineStyle = '--';
        rest = strrep(spec, '--', '');
    elseif ~isempty(strfind(spec, '-.'))
        lineStyle = '-.';
        rest = strrep(spec, '-.', '');
    elseif ~isempty(strfind(spec, ':'))
        lineStyle = ':';
        rest = strrep(spec, ':', '');
    elseif ~isempty(strfind(spec, '-'))
        lineStyle = '-';
        rest = strrep(spec, '-', '');
    else
        lineStyle = '-';
        rest = spec;
    end

    color = 'k';
    marker = 'none';
    for j = 1:length(rest)
        ch = rest(j);
        if ~isempty(strfind(colors, ch))
            color = ch;
        else
            for mk = 1:numel(markers)
                if strcmp(ch, markers{mk})
                    marker = ch;
                    break;
                end
            end
        end
    end

    st = struct('LineStyle', lineStyle, 'Marker', marker, 'Color', color);
end
