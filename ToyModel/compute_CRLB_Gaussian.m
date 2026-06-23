function CRLB_pos = compute_CRLB_Gaussian(F, H, Q, R_base, N, P0)
% 基于高斯假设（忽略离群点）计算位置 RMSE 的 CRLB 下界
% 输入：
%   F, H, Q - 系统矩阵
%   R_base  - 基础观测噪声协方差（标量或方阵）
%   N       - 时间步数
%   P0      - 初始协方差
% 输出：
%   CRLB_pos - N×1 向量，每个时刻的位置 RMSE 下界

m = size(H,1);
if isscalar(R_base)
    R = R_base * eye(m);
else
    R = R_base;
end
n = size(F,1);
P = P0;
CRLB_pos = zeros(N,1);

for k = 1:N
    P_pred = F * P * F' + Q;
    S = H * P_pred * H' + R;
    K = P_pred * H' / S;
    P = (eye(n) - K * H) * P_pred;
    P_pos = P(1:2, 1:2);
    CRLB_pos(k) = sqrt(trace(P_pos));
end
end