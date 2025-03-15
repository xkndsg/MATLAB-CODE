% function X = soft_thresholding(X,tau)
%     [U, S, V] = svd(X, 'econ') ;
%     S = sign(S) .* max(abs(S) - tau, 0) ;
%     X = U*S*V' ;
% end


function out = soft_thresholding(X, tau, svp)
    [m, n] = size(X);
    d = min(m, n); % d = min{s, n}

    % 初始奇异值计算数量
    if nargin < 3
        svp = 100;
    end
    % svp = rank(M);
    % 使用 PROPACK 或者 svds 计算前 svp 个奇异值
    opts = struct('MaxIterations', 500, 'Tolerance', 1e-6);
    [U, S, V] = svds(X, svp, 'largest', opts);

    % 软阈值化
    S = diag(S);
    S_threshold = sign(S) .* max(abs(S) - tau, 0);

    % 保留非零的奇异值
    non_zero_count = sum(S_threshold > 0);

    % 动态调整下次迭代的奇异值计算数量
    if non_zero_count < svp
        out.svp = min(non_zero_count + 1, d);
    else
        out.svp = min(non_zero_count + round(0.04 * d), d);
    end

    % 重构矩阵
    out.X = U(:, 1:non_zero_count) * diag(S_threshold(1:non_zero_count)) * V(:, 1:non_zero_count)';
end
