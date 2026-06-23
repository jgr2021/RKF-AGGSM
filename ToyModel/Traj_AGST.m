function [z, pos] = Traj_AGST(x0, F, H, Q, R_g, nu_t, sigma_t, N)
% 生成状态轨迹和观测（观测噪声 = 高斯 + 学生 t 独立相加）
% 输入：
%   x0        - 初始状态 (n×1)
%   F, H, Q   - 系统矩阵
%   R_g       - 高斯噪声协方差标量（各向同性）
%   nu_t      - 学生 t 自由度
%   sigma_t   - 学生 t 尺度标量（各向同性）
%   N         - 步数
% 输出：
%   z   - N×2 观测序列
%   pos - N×n 真实状态序列

    n = length(x0);
    m = size(H,1);
    z = zeros(N, m);
    pos = zeros(N, n);
    x = x0;
    
    for i = 1:N
        x = F * x + mvnrnd(zeros(n,1), Q)';
        pos(i,:) = x';
        
        % 高斯噪声
        noise_g = mvnrnd(zeros(1,m), R_g * eye(m))';
        % 学生 t 噪声（利用 t 分布与尺度矩阵的关系：cov = sigma_t * nu_t/(nu_t-2) 当 nu_t>2）
        noise_t = sqrt(sigma_t) * trnd(nu_t, m, 1);
        
        z(i,:) = (H * x + noise_g + noise_t)';
    end
end