function loss = test_demoRKFGSTM_general(z, pos, R_nominal, q, dt, g0, e0, w, v)
% RKF-GSTM: 带有伯努利-伽马层次先验的变分贝叶斯鲁棒卡尔曼滤波
% 输入:
%   z        - N×2 观测序列
%   pos      - N×4 真实状态（位置+速度）
%   R_nominal- 标量，名义观测噪声方差（各向同性）
%   q        - 过程噪声强度系数
%   dt       - 采样间隔
%   g0, e0   - Beta 先验形状参数
%   w, v     - Gamma 先验尺度参数

    m = 2; n = 4;
    N = size(z,1);
    
    % 系统矩阵
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
    
    % 初始状态
    x0 = [0;0;0;0];
    Xa = x0;
    Pa = eye(4);
    
    Xs = zeros(N,4);
    Ps = zeros(4,4,N);
    
    % 超参数
    % g0, e0, w, v 由外部传入
    
    for i = 1:N
        % 预测
        Xp = F * Xa;
        Pp = F * Pa * F' + Q;
        
        % 变分迭代初始化
        E_xi = 1; E_lambda = 1;
        E_log_xi = 0; E_log_lambda = 0;
        E_t = 1; E_y = 1;
        E_log_tau = psi(g0) - psi(1);
        E_log_tau_ = psi(1-g0) - psi(1);
        E_log_pi = psi(e0) - psi(1);
        E_log_pi_ = psi(1-e0) - psi(1);
        
        for j = 1:5   % 迭代次数
            % ----- 数值稳定性防御 -----
            E_t_safe = max(min(E_t, 1-1e-10), 1e-10);
            E_y_safe = max(min(E_y, 1-1e-10), 1e-10);
            E_xi_safe = max(E_xi, 1e-10);
            E_lambda_safe = max(E_lambda, 1e-10);
            
            denom_P = E_t_safe + (1 - E_t_safe) * E_xi_safe;
            denom_R = E_y_safe + (1 - E_y_safe) * E_lambda_safe;
            if denom_P < 1e-10, denom_P = 1e-10; end
            if denom_R < 1e-10, denom_R = 1e-10; end
            
            Pp_modified = (1 / denom_P) * Pp;
            R_modified = (1 / denom_R) * R;
            
            % 确保 S 矩阵对称正定
            S = H * Pp_modified * H' + R_modified;
            S = (S + S') / 2;
            [V, D] = eig(S);
            d = diag(D);
            d(d < 1e-10) = 1e-10;
            S = V * diag(d) * V';
            
            % 卡尔曼更新
            K = Pp_modified * H' / S;
            Xa = Xp + K * (z(i,:)' - H * Xp);
            Pa = (eye(n) - K * H) * Pp_modified;
            
            % 计算辅助矩阵
            A = Pa + (Xa - Xp)*(Xa - Xp)';
            B = H*Pa*H' + (z(i,:)' - H*Xa)*(z(i,:)' - H*Xa)';
            
            % 更新 t 的后验概率（使用 log-sum-exp 技巧防溢出）
            log_Pr_t1 = E_log_tau - 0.5 * trace(A / Pp);
            log_Pr_t0 = E_log_tau_ + 0.5*n*E_log_xi - 0.5*E_xi * trace(A / Pp);
            max_log_t = max(log_Pr_t1, log_Pr_t0);
            Pr_t1 = exp(log_Pr_t1 - max_log_t);
            Pr_t0 = exp(log_Pr_t0 - max_log_t);
            sum_t = Pr_t1 + Pr_t0;
            if sum_t < 1e-300
                Pr_t1 = 0.5; Pr_t0 = 0.5;
            else
                Pr_t1 = Pr_t1 / sum_t; Pr_t0 = Pr_t0 / sum_t;
            end
            
            % 更新 y 的后验概率
            log_Pr_y1 = E_log_pi - 0.5 * trace(B / R);
            log_Pr_y0 = E_log_pi_ + 0.5*m*E_log_lambda - 0.5*E_lambda * trace(B / R);
            max_log_y = max(log_Pr_y1, log_Pr_y0);
            Pr_y1 = exp(log_Pr_y1 - max_log_y);
            Pr_y0 = exp(log_Pr_y0 - max_log_y);
            sum_y = Pr_y1 + Pr_y0;
            if sum_y < 1e-300
                Pr_y1 = 0.5; Pr_y0 = 0.5;
            else
                Pr_y1 = Pr_y1 / sum_y; Pr_y0 = Pr_y0 / sum_y;
            end
            
            E_t = Pr_t1;
            E_y = Pr_y1;
            
            % 更新 xi 和 lambda（伽马分布参数）
            eta = 0.5*n*(1 - E_t) + 0.5*w;
            theta_val = 0.5*trace(A / Pp)*(1 - E_t) + 0.5*w;
            theta_val = max(theta_val, 1e-10);   % 防止零或负
            
            alpha = 0.5*m*(1 - E_y) + 0.5*v;
            beta_val = 0.5*trace(B / R)*(1 - E_y) + 0.5*v;
            beta_val = max(beta_val, 1e-10);
            
            E_xi = eta / theta_val;
            E_lambda = alpha / beta_val;
            E_log_xi = psi(max(eta,1e-10)) - log(theta_val);
            E_log_lambda = psi(max(alpha,1e-10)) - log(beta_val);
            
            % 更新 tau 和 pi（Beta 分布参数）
            g = max(g0 + E_t, 1e-10);
            h = max(2 - g0 - E_t, 1e-10);
            e = max(e0 + E_y, 1e-10);
            f = max(2 - e0 - E_y, 1e-10);
            
            E_log_tau = psi(g) - psi(g+h);
            E_log_tau_ = psi(h) - psi(g+h);
            E_log_pi = psi(e) - psi(e+f);
            E_log_pi_ = psi(f) - psi(e+f);
        end
        
        Xs(i,:) = Xa';
        Ps(:,:,i) = Pa;
    end
    
    loss = calculate_loss(pos, Xs, Ps);
end
% function loss = test_demoRKFGSTM_general(z, pos, R_nominal, q, dt, g0, e0, w, v)
% % RKF-GSTM: 带有伯努利-伽马层次先验的变分贝叶斯鲁棒卡尔曼滤波
% % 输入:
% %   z        - N×2 观测序列
% %   pos      - N×4 真实状态（位置+速度）
% %   R_nominal- 标量，名义观测噪声方差（各向同性）
% %   q        - 过程噪声强度系数
% %   dt       - 采样间隔
% 
%     m = 2; n = 4;
%     N = size(z,1);
% 
%     % 系统矩阵
%     F = [1, 0, dt, 0;
%          0, 1, 0, dt;
%          0, 0, 1, 0;
%          0, 0, 0, 1];
%     H = [1, 0, 0, 0;
%          0, 1, 0, 0];
%     Q = q * [dt^3/3, 0, dt^2/2, 0;
%              0, dt^3/3, 0, dt^2/2;
%              dt^2/2, 0, dt, 0;
%              0, dt^2/2, 0, dt];
%     R = R_nominal * eye(m);
% 
%     % 初始状态
%     x0 = [0;0;0;0];
%     Xa = x0;
%     Pa = eye(4);
% 
%     Xs = zeros(N,4);
%     Ps = zeros(4,4,N);
% 
%     % 超参数（与原代码一致）
%     % g0 = 0.85; e0 = 0.85;  % Beta 先验形状参数
%     % w = 5; v = 5;          % Gamma 先验尺度参数
% 
%     for i = 1:N
%         % 预测
%         Xp = F * Xa;
%         Pp = F * Pa * F' + Q;
% 
%         % 变分迭代初始化
%         E_xi = 1; E_lambda = 1;
%         E_log_xi = 0; E_log_lambda = 0;
%         E_t = 1; E_y = 1;
%         E_log_tau = psi(g0) - psi(1);
%         E_log_tau_ = psi(1-g0) - psi(1);
%         E_log_pi = psi(e0) - psi(1);
%         E_log_pi_ = psi(1-e0) - psi(1);
% 
%         for j = 1:5   % 迭代次数
%             % 修正协方差
%             Pp_modified = 1 / (E_t + (1-E_t)*E_xi) * Pp;
%             R_modified = 1 / (E_y + (1-E_y)*E_lambda) * R;
% 
%             % 卡尔曼更新
%             K = Pp_modified * H' / (H * Pp_modified * H' + R_modified);
%             Xa = Xp + K * (z(i,:)' - H * Xp);
%             Pa = (eye(n) - K * H) * Pp_modified;
% 
%             % 计算辅助矩阵
%             A = Pa + (Xa - Xp)*(Xa - Xp)';
%             B = H*Pa*H' + (z(i,:)' - H*Xa)*(z(i,:)' - H*Xa)';
% 
%             % 更新 t 和 y 的后验概率（伯努利变量）
%             log_Pr_t1 = E_log_tau - 0.5 * trace(A / Pp);
%             log_Pr_t0 = E_log_tau_ + 0.5*n*E_log_xi - 0.5*E_xi * trace(A / Pp);
%             Pr_t1 = exp(log_Pr_t1); Pr_t0 = exp(log_Pr_t0);
%             sum_t = Pr_t1 + Pr_t0;
%             Pr_t1 = Pr_t1 / sum_t; Pr_t0 = Pr_t0 / sum_t;
% 
%             log_Pr_y1 = E_log_pi - 0.5 * trace(B / R);
%             log_Pr_y0 = E_log_pi_ + 0.5*m*E_log_lambda - 0.5*E_lambda * trace(B / R);
%             Pr_y1 = exp(log_Pr_y1); Pr_y0 = exp(log_Pr_y0);
%             sum_y = Pr_y1 + Pr_y0;
%             Pr_y1 = Pr_y1 / sum_y; Pr_y0 = Pr_y0 / sum_y;
% 
%             E_t = Pr_t1;
%             E_y = Pr_y1;
% 
%             % 更新 xi 和 lambda（伽马分布）
%             eta = 0.5*n*(1 - E_t) + 0.5*w;
%             theta_val = 0.5*trace(A / Pp)*(1 - E_t) + 0.5*w;
%             alpha = 0.5*m*(1 - E_y) + 0.5*v;
%             beta_val = 0.5*trace(B / R)*(1 - E_y) + 0.5*v;
% 
%             E_xi = eta / theta_val;
%             E_lambda = alpha / beta_val;
%             E_log_xi = psi(eta) - log(theta_val);
%             E_log_lambda = psi(alpha) - log(beta_val);
% 
%             % 更新 tau 和 pi（Beta 分布）
%             g = g0 + E_t;
%             h = 2 - g0 - E_t;
%             e = e0 + E_y;
%             f = 2 - e0 - E_y;
% 
%             E_log_tau = psi(g) - psi(g+h);
%             E_log_tau_ = psi(h) - psi(g+h);
%             E_log_pi = psi(e) - psi(e+f);
%             E_log_pi_ = psi(f) - psi(e+f);
%         end
% 
%         Xs(i,:) = Xa';
%         Ps(:,:,i) = Pa;
%     end
% 
%     loss = calculate_loss(pos, Xs, Ps);
% end