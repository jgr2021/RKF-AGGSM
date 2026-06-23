function loss = RKF_TAU_UrbanNav(z, pos, u, dt, R, tau, c)
% RKF‑TAU 滤波器 (2D位置，带控制输入)
% z: N×2 观测, pos: N×2 真实位置, u: N×2 控制输入, dt: 采样间隔
% R: 2×2 观测噪声协方差, tau: 散度参数, c: 鲁棒控制常数

    [N, m] = size(z);
    n = 2;
    F = eye(2);
    H = eye(2);
    Q = 1e-5 * eye(2);   % 过程噪声
    Pa = 100 * eye(2);    % 初始协方差
    Xa = z(1,:)';         % 初始状态

    B = sqrtm(Q);         % 2×2
    D = sqrtm(R);         % 2×2

    Xs = zeros(N, n);
    Ps = zeros(n, n, N);

    for i = 2:N
        % 预测
        Xp = F * Xa + u(i-1,:)';
        Pp = F * Pa * F' + Q;

        % 鲁棒增益 (针对位置块)
        G = (Pp + B*D')' / (Pp + D*D');

        % 更新
        innovation = z(i,:)' - H * Xp;
        Xa = Xp + G * innovation;

        % 后验协方差 (先按标准 KF 计算，再进行收缩)
        Pa = (eye(n) - G * H) * Pp;

        % 鲁棒收缩
        Pp_sym = (Pp + Pp') / 2;   % 强制对称
        theta_val = solve_theta(Pp_sym, c, tau);
        Pa = update_V(Pp_sym, theta_val, tau);

        Xs(i,:) = Xa';
        Ps(:,:,i) = Pa;
    end
    loss = calculate_loss(pos, Xs, Ps);
end

% ------------------------------------------------------------------------
function theta = solve_theta(P, c, tau)
    lambda_max = max(eig(P));
    if tau < 1
        theta_max_valid = 0.99 / ((1-tau)*lambda_max + eps);
    else
        theta_max_valid = 1e3;
    end

    g0 = gamma_tau(P, 0, tau);
    g_max = gamma_tau(P, theta_max_valid, tau);

    if c <= g0
        theta = 0;
        return;
    elseif c >= g_max
        theta = theta_max_valid;
        return;
    end

    func = @(th) gamma_tau(P, th, tau) - c;
    options = optimset('TolX', 1e-12, 'Display', 'off');
    theta = fzero(func, [eps, theta_max_valid], options);
end

% ------------------------------------------------------------------------
function val = gamma_tau(P, theta, tau)
    lambda = eig(P);
    lambda = max(lambda(:), 0);   % 确保非负
    n = length(lambda);

    if tau == 0
        if any(theta * lambda >= 1)
            val = inf;
            return;
        end
        val = sum(log(1 - theta*lambda)) + sum(1./(1 - theta*lambda)) - n;
    elseif tau == 1
        val = sum(exp(theta*lambda) .* (theta*lambda - 1)) + n;
    else
        alpha = theta * (1 - tau);
        if any(alpha * lambda >= 1)
            val = inf;
            return;
        end
        term1 = (1 - alpha*lambda).^(tau/(tau-1));
        term2 = (1 - alpha*lambda).^(1/(tau-1));
        val = sum( -1/(tau*(1-tau))*term1 + 1/(1-tau)*term2 ) + n/tau;
    end
end

% ------------------------------------------------------------------------
function V = update_V(P, theta, tau)
    P = (P + P')/2;   % 强制对称
    try
        L = chol(P, 'lower');
    catch
        % 正定修复
        [U, D] = eig(P);
        d = diag(D);
        d(d < 1e-10) = 1e-10;
        P_fixed = U * diag(d) * U';
        P_fixed = (P_fixed + P_fixed')/2;
        L = chol(P_fixed, 'lower');
        warning('协方差矩阵不正定，已自动修复');
    end

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