function [out1, out2, I_sampled_rgb, snr1, snr2, ssim1, ssim2] = recover_image(I, is_gray)
    % 转换为双精度
    I_double = double(I);
    rng(1);
    
    % 设定采样率
    p = 0.3;
    
    % 获取图像尺寸
    [m, n, c] = size(I_double);
    
    % 生成采样掩码
    index = false(m, n);
    numTrue = round(p * m * n);
    index(randperm(m * n, numTrue)) = true;
    
    if is_gray
        % 处理灰度图像
        I_sampled = I_double .* index;
        out1 = DRSM(I_sampled, index, 0.1, 0.3, 1.6, sqrt(m*n)*p*(1-p), 2000, 1e-2);
        out2 = ADMM(I_sampled, index, sqrt(m*n)*p*(1-p), 0.6, 2000, 1e-2);
        I_sampled_rgb = I_sampled;  % 采样后的灰度图像
    else
        % 处理彩色图像
        index_re = repmat(index, 1, c);
        I_sampled = reshape(I_double, m, n * c) .* index_re;
     
        out1 = DRSM(I_sampled, index_re, 0.1, 0.3, 1.6, sqrt(m*n)*p*(1-p), 2000, 1e-5);
        out2 = ADMM(I_sampled, index_re, sqrt(m*n)*p*(1-p), 0.6, 2000, 1e-5);
        
        % 构造采样后的彩色图像
        I_sampled_rgb = reshape(I_sampled, m, n, c);

    end
    
    % 计算 SNR 和 SSIM
    snr1 = compute_snr(I_double, reshape(out1.X, size(I_double)));
    snr2 = compute_snr(I_double, reshape(out2.X, size(I_double)));
    ssim1 = compute_ssim(I_double, reshape(out1.X, size(I_double)));
    ssim2 = compute_ssim(I_double, reshape(out2.X, size(I_double)));
end

