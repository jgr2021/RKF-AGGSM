% generate_GNSS_Traj_figure.m
% 生成 UrbanNav 轨迹对比图 (GNSS_Traj.eps)
clear; clc; close all;

% 数据文件夹路径 (修改为您自己的路径)
data_dir = fullfile(pwd, 'data');

% 加载数据
tmp = load(fullfile(data_dir, 'pos.mat'));
pos = tmp.pos;
tmp = load(fullfile(data_dir, 'traj.mat'));
traj = tmp.traj;

% 如果坐标是三维的（经度、纬度、高度），此处我们已经之前处理过只取前两列，
% 但在原始 UrbanNav 数据中，您可能是用 lla2enu 转换后的平面坐标。
% 假设 pos 和 traj 已经是 N×2 平面坐标。

% 绘图
figure('Position', [100, 100, 800, 600]);
plot(traj(:,1), traj(:,2), 'r.', 'MarkerSize', 3);   % GNSS 观测点
hold on;
plot(pos(:,1), pos(:,2), 'b-', 'LineWidth', 1.5);     % 地面真值
xlabel('East (m)', 'FontSize', 12);
ylabel('North (m)', 'FontSize', 12);
legend('GNSS Measurements', 'Ground Truth', 'Location', 'best');
title('UrbanNav Trajectory');
grid on; axis equal;

% 保存为 EPS
saveas(gcf, 'GNSS_Traj.eps', 'epsc');
disp('GNSS_Traj.eps 已生成');