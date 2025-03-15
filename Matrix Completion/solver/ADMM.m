function out = ADMM(M, index, tau,rho, max_iter, tol)
    % admm_matrix_completion: 使用ADMM算法进行矩阵填充
    %
    % 输入参数:
    %   M - 输入的观测矩阵，使用 NaN 表示缺失值
    %   index - 与 M 大小相同的二值掩码矩阵，指示观测到的条目
    %   rho - ADMM 的惩罚参数
    %   max_iter - 最大迭代次数
    %   tol - 收敛的容差
    %
    % 输出:
    %   X- 补全后的矩阵

    % 初始化变量
    [m, n] = size(M);
    X = zeros(m, n);  % 初始化X
    Z = zeros(m, n);  % 初始化Z
    Y = zeros(m, n);  % 初始化Y
    out.time = 0;
    sv = 100;
    % 迭代更新
    for k = 1:max_iter
        % 保存前一次迭代的Z值
        Z_pre = Z;
        X_pre = X;
        % 更新X
        tic;
        X = project_to_observed(Z - Y/rho, M, index);

        % 更新Z：奇异值软阈值化
        x=X + Y/rho;
        % Z = soft_thresholding(x, 1/rho,100+r+m/100);
        
        out1 = soft_thresholding(x, tau/rho,sv);
        Z = out1.X;
        sv = out1.svp;
        % 更新Y：拉格朗日乘子
        Y = Y + rho * (X - Z);
        
        % 计算收敛条件
        e1 = norm(X_pre - X, 'fro')/norm(X,"fro");
        e2 = norm(Z_pre - Z, 'fro')/norm(Z,"fro");
        error = max(e1,e2);
        out.err(k) = error;
        out.X = X;
        t=toc;
        out.time = out.time+t;
        if error < tol
           out.iter = k;
            break;
        end
    end
    % out.X = X;
end


