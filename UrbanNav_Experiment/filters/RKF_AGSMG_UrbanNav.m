function loss = RKF_AGSMG_UrbanNav(z, pos, u, dt, n_Q, n_R, n_dof)
% RKF‑AGGSM for 2D position model
% n_R -> scale R matrix, n_Q -> offset coefficient, n_dof -> DOF
    [N, m] = size(z);
    n = 2; F = eye(2); H = eye(2); G = eye(2);
    Q = 1e-5 * eye(2);
    Pa = 100 * eye(2);
    Xa = z(1,:)';

    R = n_R * eye(2);
    offset = n_Q / n_R;   % 计算 d (但你的代码中 offset = n_Q/n_R)
    alpha = n_dof / 2;
    beta = n_dof / 2;
    u_ = 5/2;  U = u_ * R;

    Xs = zeros(N, n); Ps = zeros(n, n, N);

    for i = 2:N
        Xp = F * Xa + G * u(i-1,:)';
        Pp = F * Pa * F' + Q;
        E_lambda = 1;
        E_R = R;  % 您的代码初始化 E_R = U/(u_-m-1) 但随后修改，这里直接用 R
        E_lambda_ = 1;

        for jj = 1:2
            R_ = (offset + E_lambda) * E_R;
            K = Pp * H' / (H * Pp * H' + R_);
            Xa = Xp + K * (z(i,:)' - H * Xp);
            Pa = (eye(n) - K * H) * Pp;

            B = (z(i,:)' - H*Xa) * (z(i,:)' - H*Xa)' + H * Pa * H';
            k = -trace(B / E_R);
            a_coef = m + 2*alpha + 2;
            b_coef = k + offset*(m + 4*alpha + 4) - 2*beta;
            c_coef = (2*alpha + 2)*offset^2 - 4*beta*offset;
            d_coef = -2*offset^2*beta;

            if offset == 0
                E_lambda = -b_coef / a_coef;
            else
                E_lambda_ = decide_3root(solve3Polynomial(a_coef,b_coef,c_coef,d_coef), E_lambda_);
                E_lambda = E_lambda_;
            end

            D = E_lambda * B;
            uk = u_ + 1;
            Uk = U + D;
            E_R = Uk / (uk - 1 - m);
        end
        Xs(i,:) = Xa'; Ps(:,:,i) = Pa;
    end
    loss = calculate_loss(pos, Xs, Ps);
end


% 将以下全部内容追加到 RKF_AGSMG_UrbanNav.m 末尾

% ------------------------------------------------------------------------
function res = decide_3root(x, v)
    res = [];
    for i = 1:3
        if isreal(x(i)) && x(i) >= 0
            res = [res, x(i)];
        end
    end
    if length(res) ~= 1
        [~, idx] = min(abs(res - v));
        res = res(idx);
    end
end

% ------------------------------------------------------------------------
function x = solve3Polynomial(a, b, c, d)
    A = b*b - 3*a*c;   if abs(A) < 1e-14; A = 0; end
    B = b*c - 9*a*d;   if abs(B) < 1e-14; B = 0; end
    C = c*c - 3*b*d;   if abs(C) < 1e-14; C = 0; end
    DET = B*B - 4*A*C; if abs(DET) < 1e-14; DET = 0; end
    
    if (A == 0) && (B == 0)
        x1 = -c/b; x2 = x1; x3 = x1;
    elseif DET > 0
        Y1 = A*b + 1.5*a*(-B + sqrt(DET));
        Y2 = A*b + 1.5*a*(-B - sqrt(DET));
        y1 = nthroot(Y1,3);  y2 = nthroot(Y2,3);
        x1 = (-b-y1-y2)/(3*a);
        vec1 = (-b + 0.5*(y1 + y2))/(3*a);
        vec2 = 0.5*sqrt(3)*(y1 - y2)/(3*a);
        x2 = complex(vec1, vec2);
        x3 = complex(vec1, -vec2);
    elseif DET == 0 && (A ~= 0) && (B ~= 0)
        K = (b*c-9*a*d)/(b*b - 3*a*c); K = round(K,14);
        x1 = -b/a + K;   x2 = -0.5*K;   x3 = x2;
    elseif DET < 0
        sqA = sqrt(A);
        T = (A*b - 1.5*a*B)/(A*sqA);
        theta = acos(T);
        csth  = cos(theta/3);
        sn3th = sqrt(3)*sin(theta/3);
        x1 = (-b - 2*sqA*csth)/(3*a);
        x2 = (-b + sqA*(csth + sn3th))/(3*a);
        x3 = (-b + sqA*(csth - sn3th))/(3*a);
    end
    x = [x1; x2; x3];
end