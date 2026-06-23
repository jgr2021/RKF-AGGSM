%% UrbanNav 滤波器超参数扫描 (Huber, GSTM, TAU)
clear; clc; close all;
addpath('filters');

data_dir = fullfile(pwd, 'data');
[z, pos, u, dt] = load_UrbanNav_data(data_dir);

% 预先拟合高斯和学生t参数（用于GSTM和Huber的参考协方差）
err = z - pos;
R_Gaussian = diag(var(err));
R_St_mean = mean(diag(R_Gaussian)); % 简化标量，用于初始化

% 搜索网格（可自行调整）
huber_list = [0.5, 1.0, 1.345, 2.0, 3.0];
g0_list = [0.5, 0.7, 0.85, 0.95];
e0_list = [0.5, 0.7, 0.85, 0.95];
w_list  = [1, 3, 5, 10];
v_list  = [1, 3, 5, 10];
tau_list = [0, 0.3, 0.5, 0.7, 1.0];
c_list   = [0.1, 0.3, 0.5, 0.7, 1.0];

fprintf('========== 开始超参数扫描 ==========\n');

% ---- 1. 调优 Huber 阈值 ----
fprintf('\n--- 扫描 Huber 阈值 ---\n');
best_huber_rmse = inf;
best_huber_th = huber_list(1);
for th = huber_list
    try
        loss = Huber_UrbanNav(z, pos, u, dt, R_Gaussian, th);
        rmse = loss(end,1);
        fprintf('  Huber th=%.3f, RMSE=%.4f\n', th, rmse);
        if rmse < best_huber_rmse
            best_huber_rmse = rmse;
            best_huber_th = th;
        end
    catch
        fprintf('  Huber th=%.3f 失败，跳过\n', th);
    end
end
fprintf('最佳 Huber 阈值: %.3f (RMSE=%.4f)\n', best_huber_th, best_huber_rmse);

% ---- 2. 调优 GSTM 超参数 (g0, e0, w, v) ----
fprintf('\n--- 扫描 GSTM 超参数 ---\n');
best_gstm_rmse = inf;
best_gstm = struct('g0',0.85,'e0',0.85,'w',5,'v',5);
total_gstm = length(g0_list)*length(e0_list)*length(w_list)*length(v_list);
count = 0;
for g0 = g0_list
    for e0 = e0_list
        for w = w_list
            for v = v_list
                count = count + 1;
                fprintf('  GSTM %d/%d: g0=%.2f, e0=%.2f, w=%.1f, v=%.1f', count, total_gstm, g0, e0, w, v);
                try
                    loss = RKF_GSTM_UrbanNav(z, pos, u, dt, R_Gaussian, g0, e0, w, v);
                    rmse = loss(end,1);
                    fprintf(' RMSE=%.4f\n', rmse);
                    if rmse < best_gstm_rmse
                        best_gstm_rmse = rmse;
                        best_gstm.g0 = g0; best_gstm.e0 = e0; best_gstm.w = w; best_gstm.v = v;
                    end
                catch
                    fprintf(' 失败\n');
                end
            end
        end
    end
end
fprintf('最佳 GSTM: g0=%.2f, e0=%.2f, w=%.1f, v=%.1f (RMSE=%.4f)\n', ...
    best_gstm.g0, best_gstm.e0, best_gstm.w, best_gstm.v, best_gstm_rmse);

% ---- 3. 调优 TAU 超参数 (tau, c) ----
fprintf('\n--- 扫描 TAU 超参数 ---\n');
best_tau_rmse = inf;
best_tau = struct('tau',0,'c',0.5);
for tau = tau_list
    for c = c_list
        fprintf('  TAU: tau=%.1f, c=%.2f', tau, c);
        try
            loss = RKF_TAU_UrbanNav(z, pos, u, dt, R_Gaussian, tau, c);
            rmse = loss(end,1);
            fprintf(' RMSE=%.4f\n', rmse);
            if rmse < best_tau_rmse
                best_tau_rmse = rmse;
                best_tau.tau = tau; best_tau.c = c;
            end
        catch
            fprintf(' 失败\n');
        end
    end
end
fprintf('最佳 TAU: tau=%.1f, c=%.2f (RMSE=%.4f)\n', ...
    best_tau.tau, best_tau.c, best_tau_rmse);

% 保存结果
save('best_hyperparameters.mat', 'best_huber_th', 'best_gstm', 'best_tau');
fprintf('超参数已保存至 best_hyperparameters.mat\n');