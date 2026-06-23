function loss = RKF_GSTM_UrbanNav(z, pos, u, dt, R_nominal, g0, e0, w_hyper, v_hyper)
    [N, m] = size(z);
    n = 2; F = eye(2); H = eye(2); G = eye(2);
    Q = 1e-5 * eye(2);
    Pa = 100 * eye(2);
    Xa = z(1,:)';
    R = R_nominal;  % 观测噪声协方差矩阵，或可直接用标量扩成 diag？您原代码中 R 是对角阵
    % 如果 R_nominal 是标量，则转换为矩阵：
    if isscalar(R_nominal)
        R = R_nominal * eye(m);
    end

    Xs = zeros(N, n); Ps = zeros(n, n, N);

    for i = 2:N
        Xp = F * Xa + G * u(i-1,:)';
        Pp = F * Pa * F' + Q;

        E_xi = 1; E_lambda = 1;
        E_log_xi = 0; E_log_lambda = 0;
        E_t = 1; E_y = 1;
        E_log_tau = psi(g0) - psi(1);
        E_log_tau_ = psi(1-g0) - psi(1);
        E_log_pi = psi(e0) - psi(1);
        E_log_pi_ = psi(1-e0) - psi(1);

        for j = 1:2
            % 防御数值稳定性
            E_t_safe = max(min(E_t, 1-1e-10), 1e-10);
            E_y_safe = max(min(E_y, 1-1e-10), 1e-10);
            E_xi_safe = max(E_xi, 1e-10);
            E_lambda_safe = max(E_lambda, 1e-10);

            denom_P = E_t_safe + (1-E_t_safe)*E_xi_safe;
            denom_R = E_y_safe + (1-E_y_safe)*E_lambda_safe;
            Pp_mod = (1 / denom_P) * Pp;
            R_mod = (1 / denom_R) * R;

            K = Pp_mod * H' / (H * Pp_mod * H' + R_mod);
            Xa = Xp + K * (z(i,:)' - H * Xp);
            Pa = (eye(n) - K * H) * Pp_mod;

            A = Pa + (Xa - Xp)*(Xa - Xp)';
            B = H*Pa*H' + (z(i,:)' - H*Xa)*(z(i,:)' - H*Xa)';

            % 更新 t
            log_Pr_t1 = E_log_tau - 0.5 * trace(A / Pp);
            log_Pr_t0 = E_log_tau_ + 0.5*n*E_log_xi - 0.5*E_xi * trace(A / Pp);
            max_t = max(log_Pr_t1, log_Pr_t0);
            Pr_t1 = exp(log_Pr_t1 - max_t); Pr_t0 = exp(log_Pr_t0 - max_t);
            sum_t = Pr_t1 + Pr_t0;
            if sum_t < 1e-300
                Pr_t1 = 0.5; Pr_t0 = 0.5;
            else
                Pr_t1 = Pr_t1 / sum_t; Pr_t0 = Pr_t0 / sum_t;
            end

            % 更新 y
            log_Pr_y1 = E_log_pi - 0.5 * trace(B / R);
            log_Pr_y0 = E_log_pi_ + 0.5*m*E_log_lambda - 0.5*E_lambda * trace(B / R);
            max_y = max(log_Pr_y1, log_Pr_y0);
            Pr_y1 = exp(log_Pr_y1 - max_y); Pr_y0 = exp(log_Pr_y0 - max_y);
            sum_y = Pr_y1 + Pr_y0;
            if sum_y < 1e-300
                Pr_y1 = 0.5; Pr_y0 = 0.5;
            else
                Pr_y1 = Pr_y1 / sum_y; Pr_y0 = Pr_y0 / sum_y;
            end

            E_t = Pr_t1; E_y = Pr_y1;

            eta = 0.5*n*(1-E_t) + 0.5*w_hyper;
            theta_val = 0.5*trace(A / Pp)*(1-E_t) + 0.5*w_hyper;
            alpha = 0.5*m*(1-E_y) + 0.5*v_hyper;
            beta_val = 0.5*trace(B / R)*(1-E_y) + 0.5*v_hyper;

            theta_val = max(theta_val, 1e-10);
            beta_val = max(beta_val, 1e-10);
            E_xi = eta / theta_val;
            E_lambda = alpha / beta_val;
            E_log_xi = psi(max(eta,1e-10)) - log(theta_val);
            E_log_lambda = psi(max(alpha,1e-10)) - log(beta_val);

            g = max(g0 + E_t, 1e-10);
            h = max(2 - g0 - E_t, 1e-10);
            e = max(e0 + E_y, 1e-10);
            f = max(2 - e0 - E_y, 1e-10);
            E_log_tau = psi(g) - psi(g+h);
            E_log_tau_ = psi(h) - psi(g+h);
            E_log_pi = psi(e) - psi(e+f);
            E_log_pi_ = psi(f) - psi(e+f);
        end
        Xs(i,:) = Xa'; Ps(:,:,i) = Pa;
    end
    loss = calculate_loss(pos, Xs, Ps);
end