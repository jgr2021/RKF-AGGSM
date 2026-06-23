function loss = KF_UrbanNav(z, pos, u, dt, R)
% 标准卡尔曼滤波器 (2D位置，带控制输入)
% R - 观测噪声协方差矩阵 (2x2)

    [N, m] = size(z);
    n = 2;
    F = eye(2);
    H = eye(2);
    G = eye(2);
    Q = 1e-5 * eye(2);   % 过程噪声协方差 (可调整)
    Pa = 100 * eye(2);
    Xa = z(1,:)';

    Xs = zeros(N, n);
    Ps = zeros(n, n, N);

    for i = 2:N
        % 预测: x = F*x + G*u
        Xp = F * Xa + G * u(i-1,:)';
        Pp = F * Pa * F' + Q;
        % 更新
        K = Pp * H' / (H * Pp * H' + R);
        Xa = Xp + K * (z(i,:)' - H * Xp);
        Pa = (eye(n) - K * H) * Pp;
        Xs(i,:) = Xa';
        Ps(:,:,i) = Pa;
    end
    loss = calculate_loss(pos, Xs, Ps);
end