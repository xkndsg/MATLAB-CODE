function out = resgray(I)
% 读取灰度图像
I_gray = double(I);  % 灰度图像转换为 double 类型
% 获取图像的尺寸
[M, N] = size(I_gray);
p = 0.3;
% 生成采样掩码
index = false(M, N);  % 创建一个 M x N 的全 false 布尔矩阵
numTrue = round(0.3 * M * N);  % 计算需要设置为 true 的元素数量
trueIndices = randperm(M * N, numTrue);  % 随机选择 numTrue 个唯一的索引
index(trueIndices) = true;  % 将选定索引对应的元素设置为 true
% 对灰度图像进行采样
I_sampled = I_gray .* index;
% 使用 DRSM 填充
out_DRSM = DRSM(I_sampled, index,0.1, 0.3 , 1.6,  sqrt(m*n)*p*(1-p) , 2000, 1e-6);
I_filled_DRSM = out_DRSM.X;  
% 使用 ADMM 填充
out_ADMM = ADMM(I_sampled, index, sqrt(m*n)*p*(1-p), 0.6 , 2000, 1e-6);
I_filled_ADMM = out_ADMM.X;  

disp(['DRSM 时间: ', num2str(out_DRSM.time)]);
disp(['ADMM 时间: ', num2str(out_ADMM.time)]);

% 计算恢复误差
error1 = abs(~index .* (I_gray - I_filled_DRSM));
out.error_DRSM = norm(error1,"fro")/norm(I_gray,"fro");
error2= abs(~index .* (I_gray - I_filled_ADMM));
out.error_ADMM = norm(error2,"fro")/norm(I_gray,"fro");
out.I_sampled = I_sampled;
out.I_filled_DRSM = I_filled_DRSM;
out.I_filled_ADMM = I_filled_ADMM;
out.iter_DRSM = out_DRSM.iter;
out.iter_ADMM = out_ADMM.iter;
end