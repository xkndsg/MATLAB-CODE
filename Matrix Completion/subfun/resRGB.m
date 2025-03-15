function out = resRGB(I)
I_double = double(I); % 转换为双精度
% 获取图像的尺寸
[M, N, C] = size(I_double);
% 生成采样掩码
index = false(M, N);  % 创建一个 M x N 的全 false 布尔矩阵
numTrue = round(0.3 * M * N);  % 计算需要设置为 true 的元素数量
trueIndices = randperm(M * N, numTrue);  % 随机选择 numTrue 个唯一的索引
index(trueIndices) = true;  % 将选定索引对应的元素设置为 true
% 将原图像展平为二维矩阵
I_re = reshape(I_double, M, N * C);  % 将图像矩阵展平成二维矩阵
index_re = repmat(index, 1, C);  % 将采样掩码复制到三个通道
I_sampled = I_re .* index_re;  % 对原图像进行采样
% 使用算法填充缺失值
out1 = DRSM(I_sampled, index_re, 0.1, 0.3 , 1.6,  sqrt(m*n)*p*(1-p) , 2000, 1e-5);
out2 = ADMM(I_sampled, index_re,  sqrt(m*n)*p*(1-p), 0.6 , 2000, 1e-5);
% 将填充后的矩阵分解回 RGB 图像
out.I_filled_DRSM = reshape(out1.X, M, N, C);
out.I_filled_ADMM = reshape(out2.X, M, N, C);
% 构造采样后的图像
DD_R = index .* I_double(:, :, 1);
DD_G = index .* I_double(:, :, 2);
DD_B = index .* I_double(:, :, 3);
out.I_sampled = cat(3, DD_R, DD_G, DD_B);

% 显示计算时间
disp(['DRSM time: ', num2str(out1.time)]);
disp(['ADMM time: ', num2str(out2.time)]);
% disp(['SVT time: ', num2str(out3.time)]);
% 计算恢复误差
error1 = abs(~index_re .* (I_re - out1.X));
out.error_DRSM= norm(error1, "fro")/norm(reshape(I_double,M,N*C),"fro");
error2 = abs(~index_re .* (I_re - out2.X));
out.error_ADMM = norm(error2, "fro")/norm(reshape(I_double,M,N*C),"fro");
out.iter_DRSM = out1.iter;
out.iter_ADMM = out2.iter;
end