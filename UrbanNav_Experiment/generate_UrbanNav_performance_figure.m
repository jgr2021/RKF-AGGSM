% generate_UrbanNav_performance_figure.m
% 生成 UrbanNav RMSE/MAE 曲线图（左右子图，y 轴上限减半，标记错开）
clear; clc; close all;

% 加载结果文件
load('results/UrbanNav_full_results.mat');  % 包含各滤波器的损失变量

% 时间轴 (采样间隔 dt=0.1 s)
dt = 0.1;
N = size(loss_KF, 1);
time_axis = (0:N-1) * dt;   % 从 0 开始，单位秒

% 定义线型和图例名（9 种滤波器）
line_styles = {'--^b', '--sk', '--dg', '--+m', '--xc', '-.b', '-.g', '-or', '-sm'};
legend_names = {'KF', 'PF', 'Huber', 'RKF-Gaussian', 'RKF-ST', ...
                'RKF-GSTM', 'RKF-TAU', 'RKF-AGST', 'RKF-AGSlash'};

% 提取 RMSE 和 MAE 数据（顺序对应）
plot_data_RMSE = {
    loss_KF(:,1), loss_PF(:,1), loss_Huber(:,1), ...
    loss_RKFG(:,1), loss_RKFT(:,1), ...
    loss_GSTM(:,1), loss_TAU(:,1), loss_RKFAGSMG(:,1), loss_Slash(:,1)
};
plot_data_MAE = {
    loss_KF(:,2), loss_PF(:,2), loss_Huber(:,2), ...
    loss_RKFG(:,2), loss_RKFT(:,2), ...
    loss_GSTM(:,2), loss_TAU(:,2), loss_RKFAGSMG(:,2), loss_Slash(:,2)
};

% 标记错开参数
base_step = 500;          % 基础步长（点数）
offsets = 0:8;            % 每条曲线的起始偏移
marker_indices = cell(9,1);
for i = 1:9
    start_idx = 1 + offsets(i);
    step = base_step + (i-1)*10;
    marker_indices{i} = start_idx:step:length(time_axis);
end

% 创建图形（左右子图，宽度增加，高度适中）
figure('Position', [100, 100, 1200, 500]);

% ----- 左子图：RMSE -----
subplot(1,2,1);
hold on;
for i = 1:9
    plot(time_axis, plot_data_RMSE{i}, line_styles{i}, ...
         'LineWidth', 1.5, 'MarkerSize', 5, 'MarkerIndices', marker_indices{i});
end
hold off;
xlabel('Time (s)', 'FontSize', 12);
ylabel('RMSE (m)', 'FontSize', 12);
legend(legend_names, 'Location', 'best', 'FontSize', 10);
box on; grid on;
set(gca, 'FontSize', 11);
title('UrbanNav - RMSE', 'FontSize', 12);
xlim([0, max(time_axis)]);
ylim([0, 25]);   % 原上限可能 50+，减半至 25（可根据数据实际最大值调整）
% 如果某些滤波器的 RMSE 超过 25，可适当提高，但 25 通常能包含大部分曲线

% ----- 右子图：MAE -----
subplot(1,2,2);
hold on;
for i = 1:9
    plot(time_axis, plot_data_MAE{i}, line_styles{i}, ...
         'LineWidth', 1.5, 'MarkerSize', 5, 'MarkerIndices', marker_indices{i});
end
hold off;
xlabel('Time (s)', 'FontSize', 12);
ylabel('MAE (m)', 'FontSize', 12);
legend(legend_names, 'Location', 'best', 'FontSize', 10);
box on; grid on;
set(gca, 'FontSize', 11);
title('UrbanNav - MAE', 'FontSize', 12);
xlim([0, max(time_axis)]);
ylim([0, 20]);   % MAE 通常略小于 RMSE，设上限 20 即可

% 保存为 EPS
saveas(gcf, 'UrbanNav_RMSE_MAE.eps', 'epsc');
disp('UrbanNav_RMSE_MAE.eps 已生成（左右布局，y轴上限减半，标记错开）');