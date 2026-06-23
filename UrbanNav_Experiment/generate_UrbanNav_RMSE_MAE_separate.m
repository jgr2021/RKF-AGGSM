% generate_UrbanNav_RMSE_MAE_separate.m
% 分别生成 UrbanNav 的 RMSE 和 MAE 图（独立 figure，无 RKF-AGSlash）
clear; clc; close all;

% 加载结果文件
load('results/UrbanNav_full_results.mat');  % 包含各滤波器的损失变量

% 时间轴 (采样间隔 dt=0.1 s)
dt = 0.1;
N = size(loss_KF, 1);
time_axis = (0:N-1) * dt;   % 从 0 开始，单位秒

% 定义线型和图例名（仅前8种滤波器，不含 Slash）
line_styles = {'--^b', '--sk', '--dg', '--+m', '--xc', '-.b', '-.g', '-or', '-sm'};
legend_names = {'KF', 'PF', 'Huber', 'RKF-Gaussian', 'RKF-ST', ...
                'RKF-GSTM', 'RKF-TAU', 'RKF-AGST', 'RKF-AGSlash'};

% 提取 RMSE 和 MAE 数据（顺序对应，只取前8个）
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

% 标记错开参数（8条曲线）
base_step = 500;          % 基础步长（点数）
offsets = 0:8;            % 每条曲线的起始偏移
marker_indices = cell(9,1);
for i = 1:9
    start_idx = 1 + offsets(i);
    step = base_step + (i-1)*10;
    marker_indices{i} = start_idx:step:length(time_axis);
end

% ================== 图1：RMSE ==================
figure('Name', 'UrbanNav RMSE', 'Position', [100, 100, 800, 500]);
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
ylim([0, 25]);   % 上限减半，可根据数据实际最大值调整
saveas(gcf, 'UrbanNav_RMSE.eps', 'epsc');
fprintf('UrbanNav_RMSE.eps 已生成\n');

% ================== 图2：MAE ==================
figure('Name', 'UrbanNav MAE', 'Position', [100, 100, 800, 500]);
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
ylim([0, 20]);   % MAE 上限适当调低
saveas(gcf, 'UrbanNav_MAE.eps', 'epsc');
fprintf('UrbanNav_MAE.eps 已生成\n');

disp('所有图已生成完毕（RMSE 和 MAE 分开，无 RKF-AGSlash）');