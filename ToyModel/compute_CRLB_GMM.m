function CRLB_pos = compute_CRLB_GMM(F, H, Q, weights, covs, N_steps, P0, n_fisher)
% 计算 GMM 观测噪声下的 CRLB（位置 RMSE 下界）
% 输入:
%   F, H, Q - 系统矩阵
%   weights - 1×K 混合权重
%   covs    - m×m×K 各分量协方差（假设均值为零）
%   N_steps - 仿真步数
%   P0      - 初始状态协方差
%   n_fisher- 用于估计 Fisher 信息的样本数（建议 1e5）
% 输出:
%   CRLB_pos - N_steps×1，每个时刻位置估计的 RMSE 下界

    m = size(H, 1);
    % 估计 Fisher 信息
    I_fish = estimate_gmm_fisher_zero_mean(weights, covs, n_fisher);
    R_eff = inv(I_fish);   % 等效观测噪声协方差
    
    % 递推 CRLB
    n = size(F, 1);
    P = P0;
    CRLB_pos = zeros(N_steps, 1);
    for k = 1:N_steps
        P_pred = F * P * F' + Q;
        S = H * P_pred * H' + R_eff;
        K = P_pred * H' / S;
        P = (eye(n) - K * H) * P_pred;
        P_pos = P(1:2, 1:2);
        CRLB_pos(k) = sqrt(trace(P_pos));
    end
end
function I_fish = estimate_gmm_fisher_zero_mean(weights, covs, N_samples)
    m = size(covs,1);
    K = length(weights);
    % 生成样本
    comp = randsample(1:K, N_samples, true, weights);
    samples = zeros(N_samples, m);
    for k = 1:K
        idx = (comp == k);
        n_k = sum(idx);
        samples(idx, :) = mvnrnd(zeros(1,m), covs(:,:,k), n_k);
    end
    scores = zeros(N_samples, m);
    for i = 1:N_samples
        v = samples(i,:)';
        log_pdf = zeros(K,1);
        grad_vals = zeros(m,K);
        for k = 1:K
            Sigma_k = covs(:,:,k);
            % 使用对数密度避免下溢
            Lk = chol(Sigma_k, 'lower');
            diff = Lk \ v;
            log_pdf(k) = log(weights(k)) - 0.5*sum(diff.^2) - sum(log(diag(Lk))) - (m/2)*log(2*pi);
            grad_vals(:,k) = - Sigma_k \ v;
        end
        max_log = max(log_pdf);
        pdf_vals_scaled = exp(log_pdf - max_log);   % 防止下溢
        total_scaled = sum(pdf_vals_scaled);
        s = grad_vals * (pdf_vals_scaled ./ total_scaled);
        scores(i,:) = s';
    end
    I_fish = (scores' * scores) / N_samples;
    I_fish = (I_fish + I_fish') / 2;
    [V,D] = eig(I_fish);
    d = diag(D);
    d(d < 1e-10) = 1e-10;
    I_fish = V * diag(d) * V';
end
% function I_fish = estimate_gmm_fisher_zero_mean(weights, covs, N_samples)
% % 零均值 GMM 的 Fisher 信息矩阵估计（蒙特卡洛）
%     m = size(covs, 1);
%     K = length(weights);
% 
%     % 生成样本
%     comp = randsample(1:K, N_samples, true, weights);
%     samples = zeros(N_samples, m);
%     for k = 1:K
%         idx = (comp == k);
%         n_k = sum(idx);
%         samples(idx, :) = mvnrnd(zeros(1, m), covs(:, :, k), n_k);
%     end
% 
%     % 逐样本计算得分
%     scores = zeros(N_samples, m);
%     for i = 1:N_samples
%         v = samples(i, :)';
%         pdf_vals = zeros(K, 1);
%         grad_vals = zeros(m, K);
%         for k = 1:K
%             Sigma_k = covs(:, :, k);
%             pdf_vals(k) = weights(k) * mvnpdf(v', zeros(1, m), Sigma_k);
%             grad_vals(:, k) = - Sigma_k \ v;
%         end
%         total_pdf = sum(pdf_vals);
%         if total_pdf < 1e-300
%             s = zeros(m, 1);
%         else
%             s = grad_vals * (pdf_vals ./ total_pdf);
%         end
%         scores(i, :) = s';
%     end
% 
%     % Fisher 信息 = 得分外积的均值
%     I_fish = (scores' * scores) / N_samples;
%     I_fish = (I_fish + I_fish') / 2;
%     % 确保正定
%     [V, D] = eig(I_fish);
%     d = diag(D);
%     d(d < 1e-10) = 1e-10;
%     I_fish = V * diag(d) * V';
% end