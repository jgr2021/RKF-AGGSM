function Samples = GMMrnd(weights, means, covs, N)
% weights: 1 x K 混合权重
% means:   K x D 各分量均值（通常为零向量）
% covs:    D x D x K 各分量协方差矩阵
% N:       样本数量
    K = length(weights);
    D = size(means, 2);
    Samples = zeros(N, D);
    comp = randsample(1:K, N, true, weights);
    for k = 1:K
        idx = (comp == k);
        n_k = sum(idx);
        Samples(idx, :) = mvnrnd(means(k, :), covs(:, :, k), n_k);
    end
end