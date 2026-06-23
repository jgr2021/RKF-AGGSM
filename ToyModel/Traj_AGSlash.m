function [z, pos] = Traj_AGSlash(x0, F, H, Q, R_g, alpha_s, scale_s, N)
% 生成状态轨迹和观测（观测噪声 = 高斯 + Slash 独立相加）
% 输入：
%   x0       - 初始状态 (n×1)
%   F, H, Q  - 系统矩阵
%   R_g      - 高斯噪声协方差标量（各向同性）
%   alpha_s  - Slash 分布形状参数（>0，越小尾部越重）
%   scale_s  - Slash 分布尺度参数（标量）
%   N        - 步数
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
        % Slash 噪声（各向同性，独立同分布）
        noise_slash = slashrnd(0, sqrt(scale_s), alpha_s, [m,1]);
        
        z(i,:) = (H * x + noise_g + noise_slash)';
    end
end