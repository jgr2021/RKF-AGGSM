function loss = RKF_Slash_UrbanNav(z, pos, u, dt, n_Q, n_R, n_dof)
% 输入参数与 RKF_AGSMG_UrbanNav 一致
% n_R : 观测噪声尺度因子
% n_Q : 偏移参数 (offset = n_Q / n_R)
% n_dof: 分布自由度
    [N, m] = size(z);
    n = 2;
    F = eye(2); H = eye(2); G = eye(2);
    Q = 1e-3 * eye(2);        % 适当增大过程噪声，提高鲁棒性
    Pa = 100 * eye(2);
    Xa = z(1,:)';

    R = n_R * eye(2);
    offset = n_Q / n_R;
    dof = n_dof;

    u_ = 5/2;              % 逆 Wishart 自由度（与测试版本一致）
    U = u_ * R;
    E_R = U / (u_ - m - 1);   % 正确的先验期望

    Xs = zeros(N, n);
    Ps = zeros(n, n, N);

    for i = 2:N
        Xp = F * Xa + G * u(i-1,:)';
        Pp = F * Pa * F' + Q;
        E_lambda = 1;

        for jj = 1:3
            R_ = (offset + E_lambda) * E_R;
            K = Pp * H' / (H * Pp * H' + R_);
            Xa = Xp + K * (z(i,:)' - H * Xp);
            Pa = (eye(n) - K * H) * Pp;

            B = (z(i,:)' - H*Xa)*(z(i,:)' - H*Xa)' + H * Pa * H';
            k = -trace(B / E_R);
            a = m + dof - 2;
            b = k + m*offset + 2*(dof-2)*offset;
            c = (dof-2)*offset^2;

            % 求解 a*E_lambda^2 + b*E_lambda + c = 0
            if abs(a) < 1e-12
                E_lambda = -c / b;
            else
                disc = b^2 - 4*a*c;
                if disc < 0
                    E_lambda = 1 + offset;
                else
                    E_lambda = (sqrt(disc) - b) / (2*a);
                end
            end
            if ~(isreal(E_lambda) && E_lambda > 0)
                E_lambda = 1 + offset;
            end
            E_lambda = max(1 + offset, E_lambda);

            D = E_lambda * B;
            uk = u_ + 1;
            Uk = U + D;
            E_R = Uk / (uk - 1 - m);
        end

        Xs(i,:) = Xa';
        Ps(:,:,i) = Pa;
    end

    loss = calculate_loss(pos, Xs, Ps);
end