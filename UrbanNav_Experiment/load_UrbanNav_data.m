function [z, pos, u, dt] = load_UrbanNav_data(data_dir)
% 从给定文件夹加载 pos.mat, traj.mat, data_vel.mat 和 time.mat (可选)
% 输出皆为 Nx2 矩阵，适应二维位置模型
% z   : 观测位置 (Nx2)
% pos : 真实位置 (Nx2)
% u   : 控制输入 (速度*dt, Nx2)
% dt  : 采样间隔

    % 1. 真实位置
    tmp = load(fullfile(data_dir, 'pos.mat'));
    pos_raw = tmp.pos;
    if size(pos_raw,2) >= 3
        pos = pos_raw(:, 1:2);   % 只保留前两列平面坐标
        warning('pos 包含超过2列，已自动截取前两列');
    else
        pos = pos_raw;
    end

    % 2. GPS 观测
    tmp = load(fullfile(data_dir, 'traj.mat'));
    z_raw = tmp.traj;
    if size(z_raw,2) >= 3
        z = z_raw(:, 1:2);
        warning('traj 包含超过2列，已自动截取前两列');
    else
        z = z_raw;
    end

    % 3. 速度数据 (控制量)
    tmp = load(fullfile(data_dir, 'data_vel.mat'));
    vel = tmp.data_vel;
    if size(vel,2) == 1
        vel = [vel, vel];   % 单列速度复制为两列（假设各向同性）
        warning('data_vel 只有一列，已复制为两列');
    elseif size(vel,2) > 2
        vel = vel(:, 1:2);
        warning('data_vel 超过两列，已截取前两列');
    end

    % 4. 时间戳 (可选)
    if exist(fullfile(data_dir, 'time.mat'), 'file')
        tmp = load(fullfile(data_dir, 'time.mat'));
        t = tmp.t;
    else
        t = (1:size(pos,1))';
        warning('未找到 time.mat，使用默认时间轴 dt=1s');
    end

    % 5. 对齐长度
    N = min([size(pos,1), size(z,1), size(vel,1), length(t)]);
    pos = pos(1:N, :);
    z   = z(1:N, :);
    vel = vel(1:N, :);
    t   = t(1:N);

    % 采样间隔
    dt=0.1;
    % 控制输入 u = 速度*dt
    u = vel * dt;

    % 最终维度检查
    assert(size(z,2)==2 && size(pos,2)==2 && size(u,2)==2, ...
        '输出变量必须为 Nx2，请检查原始数据');
    fprintf('数据加载完毕: %d 个样本, dt=%.2f s\n', N, dt);
end