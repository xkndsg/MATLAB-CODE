clc;
clear;
addpath("figure\","solver\","subfun\")

% 读取多个图像
image_files = {'cloth.png', 'building.png', 'window.png', 'Barbara.png', 'Cameraman.png', 'peppers.png'};
num_images = numel(image_files);

% 初始化存储变量
X_filled1_all = cell(1, num_images);
X_filled2_all = cell(1, num_images);
I_sampled_all = cell(1, num_images);
snr1_all = zeros(1, num_images);
snr2_all = zeros(1, num_images);
ssim1_all = zeros(1, num_images);
ssim2_all = zeros(1, num_images);
time1_all = zeros(1, num_images);
time2_all = zeros(1, num_images);
iter1_all = zeros(1, num_images);
iter2_all = zeros(1, num_images);

% 逐个处理图像
for i = 1:1
    I = imread(image_files{i});
    is_gray = size(I, 3) == 1;
    
    % 直接从 recover_image 获取输出
    [out1, out2, I_sampled, snr1, snr2, ssim1, ssim2] = recover_image(I, is_gray);
    
    % 存储结果
    I_sampled_all{i} = I_sampled;
    snr1_all(i) = snr1;
    snr2_all(i) = snr2;
    ssim1_all(i) = ssim1;
    ssim2_all(i) = ssim2;
    time1_all(i) = out1.time;
    time2_all(i) = out2.time;
    iter1_all(i) = out1.iter;
    iter2_all(i) = out2.iter;

    % 显示 SNR、SSIM、计算时间和迭代次数
    disp(['Image ', num2str(i), ' - DRSM: SNR = ', num2str(snr1), ' dB, SSIM = ', num2str(ssim1), ...
          ', Time = ', num2str(out1.time), 's, Iterations = ', num2str(out1.iter)]);
    disp(['Image ', num2str(i), ' - ADMM: SNR = ', num2str(snr2), ' dB, SSIM = ', num2str(ssim2), ...
          ', Time = ', num2str(out2.time), 's, Iterations = ', num2str(out2.iter)]);
end
%%
% 画图：每张图像四列（原图、采样图、DRSM 恢复、ADMM 恢复）
figure;
for i = 1:num_images
    % 原始图像
    subplot(num_images, 4, (i - 1) * 4 + 1);
    imshow(uint8(imread(image_files{i})));
    % title(['Original ', num2str(i)]);
    
    % 采样图像
    subplot(num_images, 4, (i - 1) * 4 + 2);
    imshow(uint8(I_sampled_all{i}));
    % title(['Sampled ', num2str(i)]);
    
    % DRSM 恢复图像
    subplot(num_images, 4, (i - 1) * 4 + 3);
    imshow(uint8(reshape(out1.X, size(I))));
    % title(['DRSM ', num2str(i)]);
    
    % ADMM 恢复图像
    subplot(num_images, 4, (i - 1) * 4 + 4);
    imshow(uint8(reshape(out2.X, size(I))));
    % title(['ADMM ', num2str(i)]);
end