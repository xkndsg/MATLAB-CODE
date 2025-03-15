load OR.mat

figure(1);
ha = tightsubplot(2, 4, [0.02 0.02], [0.02 0.02],[0.02 0.02]);
% 显示原始图像
axes(ha(1));
imshow(uint8(I));  

% 显示采样后的图像
axes(ha(2));
imshow(uint8(I_sampled_rgb));  

% 显示 DRSM 填充后的图像
axes(ha(3));
imshow(uint8(X_filled1));  

% 显示 ADMM 填充后的图像
axes(ha(4));
imshow(uint8(X_filled2));  

load LOW.mat
% 显示原始图像
axes(ha(5));
imshow(uint8(I_low_rank));  

% 显示采样后的图像
axes(ha(6));
imshow(uint8(I_sampled_rgb));  

% 显示 DRSM 填充后的图像
axes(ha(7));
imshow(uint8(X_filled1));  

% 显示 ADMM 填充后的图像
axes(ha(8));
imshow(uint8(X_filled2));  