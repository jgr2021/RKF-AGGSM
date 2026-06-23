function loss = test_demoRKFtau_general(z, pos, R_nominal, q, dt, tau, c)
% RKF-TAU: 基于 γ-τ 散度的鲁棒卡尔曼滤波 (适配 n=4, m=2)
% 输入:
%   z        - N×2 观测序列
%   pos      - N×4 真实状态
%   R_nominal- 标量，名义观测噪声方差
%   q        - 过程噪声强度系数
%   dt       - 采样间隔
%   tau      - 散度参数 (0 ≤ τ ≤ 1)
%   c        - 鲁棒控制常数

    m = 2; n = 4;
    N = size(z,1);
    
    F = [1, 0, dt, 0;
         0, 1, 0, dt;
         0, 0, 1, 0;
         0, 0, 0, 1];
    H = [1, 0, 0, 0;
         0, 1, 0, 0];
    Q = q * [dt^3/3, 0, dt^2/2, 0;
             0, dt^3/3, 0, dt^2/2;
             dt^2/2, 0, dt, 0;
             0, dt^2/2, 0, dt];
    R = R_nominal * eye(m);
    
    % 过程噪声的平方根：提取位置相关部分 (2×2)，因为原算法假设状态维=2
    % 我们取 Q 中与位置对应的子矩阵的 Cholesky 因子
    Q_pos = Q(1:2, 1:2);          % 位置的过程噪声协方差
    B = sqrtm(Q_pos);             % 2×2
    
    % 观测噪声平方根
    D = sqrtm(R);                 % 2×2
    
    x0 = [0;0;0;0];
    Xa = x0;
    Pa = eye(4);
    
    Xs = zeros(N,4);
    Ps = zeros(4,4,N);
    
    for i = 1:N
        % 预测
        Xp = F * Xa;
        Pp = F * Pa * F' + Q;
        
        % 提取位置相关的协方差块
        Pp_pos = Pp(1:2, 1:2);
        Pxz = Pp(1:2, 1:2);       % 观测与状态的互协方差 (H*Pp*H' 就是 Pp_pos)
        
        % 鲁棒卡尔曼增益（参考原代码形式）
        % 原代码: G = (F*Pa*F + B*D')' / (H*Pa*H' + D*D')
        % 这里 F*Pa*F' 对应预测协方差的位置块，H*Pa*H' 也是位置块
        G = (Pp_pos + B*D')' / (Pp_pos + D*D');
        % 注意：G 是 2×2
        
        % 状态更新（仅位置分量）
        innovation = z(i,:)' - H * Xp;
        Xa = Xp;
        Xa(1:2) = Xp(1:2) + G * innovation;
        
        % 后验协方差（标准卡尔曼形式，但增益 G 只作用于位置）
        K_full = zeros(n, m);
        K_full(1:2, :) = G;
        Pa = (eye(n) - K_full * H) * Pp;
        
        Pp_pos = (Pp_pos + Pp_pos')/2;  % 强制对称
        % 鲁棒协方差修正（对位置块应用 γ-τ 收缩）
        theta = solve_theta(Pp_pos, c, tau);
        Pa_pos = update_V(Pp_pos, theta, tau);
        Pa(1:2, 1:2) = Pa_pos;
        
        Xs(i,:) = Xa';
        Ps(:,:,i) = Pa;
    end
    
    loss = calculate_loss(pos, Xs, Ps);
end

% ========== 以下辅助函数保持不变 ==========
function V = update_V(P, theta, tau)
    L = chol(P, 'lower');
    n = size(P,1);
    if tau == 0
        M = eye(n) - theta * (L' * L);
        V = L * (M \ L');
    elseif tau == 1
        V = L * expm(theta * (L' * L)) * L';
    else
        M = L' * L;
        [U, D] = eig(M);
        lambda = diag(D);
        lambda = max(lambda, 0);
        alpha = theta * (1 - tau);
        beta = (1 - alpha * lambda).^(1/(tau-1));
        M_pow = U * diag(beta) * U';
        V = L * M_pow * L';
    end
end

function theta = solve_theta(P, c, tau)
% 鲁棒求解 theta 使得 gamma_tau(P, theta, tau) = c
% 若 c 大于 gamma_tau 的最大可能值，则返回定义域上界（最大鲁棒收缩）
    lambda_max = max(eig(P));
    lambda_max = max(lambda_max, 1e-6);  % 防止零特征值
    
    % 定义域上界 (确保 theta 在有效范围内)
    if tau < 1
        theta_max_valid = 0.999 / ((1 - tau) * lambda_max); % 留一点裕量
    else
        theta_max_valid = 1e3;  % tau=1 时无理论上限，取较大值
    end
    
    % 计算 gamma_tau 在 theta=0 和 theta=theta_max_valid 处的值
    g0 = gamma_tau(P, 0, tau);  % 理论上应为 0 (当 theta=0 时)
    g_max = gamma_tau(P, theta_max_valid, tau);
    
    % 如果 c 小于等于 0，直接返回 0
    if c <= g0
        theta = 0;
        return;
    end
    
    % 如果 c 大于最大值，返回最大有效 theta
    if c >= g_max
        theta = theta_max_valid;
        return;
    end
    
    % 否则在 [eps, theta_max_valid] 内搜索零点
    theta_low = eps;
    theta_high = theta_max_valid;
    
    % 使用 fzero 求解
    func = @(th) gamma_tau(P, th, tau) - c;
    options = optimset('TolX', 1e-12, 'Display', 'off');
    try
        theta = fzero(func, [theta_low, theta_high], options);
    catch
        % 如果 fzero 失败，使用网格搜索近似
        warning('fzero 失败，改用网格搜索');
        th_grid = linspace(theta_low, theta_high, 1000);
        g_vals = arrayfun(@(th) gamma_tau(P, th, tau), th_grid);
        [~, idx] = min(abs(g_vals - c));
        theta = th_grid(idx);
    end
end
function val = gamma_tau(P, theta, tau)
    lambda = eig(P);
    lambda = lambda(:);
    lambda = max(lambda, 0);  % 确保非负
    n = length(lambda);
    
    if tau == 0
        % 需要 1 - theta*lambda > 0
        valid = (theta * lambda) < 1;
        if ~all(valid)
            val = inf;  % 无效区域返回大值，促使求解器避开
            return;
        end
        val = sum(log(1 - theta*lambda)) + sum(1./(1 - theta*lambda)) - n;
        
    elseif tau == 1
        exp_terms = exp(theta * lambda);
        val = sum(exp_terms .* (theta * lambda - 1)) + n;
        
    else
        alpha = theta * (1 - tau);
        % 要求 1 - alpha*lambda > 0
        valid = (alpha * lambda) < 1;
        if ~all(valid)
            val = inf;
            return;
        end
        term1 = (1 - alpha * lambda).^(tau/(tau-1));
        term2 = (1 - alpha * lambda).^(1/(tau-1));
        val = sum( -1/(tau*(1-tau)) * term1 + 1/(1-tau) * term2 ) + n/tau;
    end
end