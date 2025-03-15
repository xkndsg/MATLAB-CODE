clc;
clear;
% 定义文件名列表
addpath("figure\","solver\","subfun\")
names = {'Barbara.bmp', 'Peppers.bmp', 'BaboonRGB.bmp'};

% 假设 resRGB 和 resgray 返回的结构体具有以下字段
out = struct('I_sampled', [], 'I_filled_DRSM', [], 'I_filled_ADMM', [], 'error_DRSM', [], 'error_ADMM', []);

% 创建紧密排列的子图布局
figure(1);
ha = tightsubplot(4, 3, [0.01 0.01], [0.01 0.01], [0.01 0.01]);

for i = 1:3
    % 读取图像
    I = imread(names{i});
    
    % 判断是否为RGB图像
    if size(I, 3) == 3  % 如果有3个通道
        out_temp = resRGB(I);  % 使用resRGB函数处理RGB图像
    else
        out_temp = resgray(I);  % 使用resgray函数处理灰度图像
    end

    % 将返回的结构体字段赋值到结构体数组中
    out(i).I_sampled = out_temp.I_sampled;
    out(i).I_filled_DRSM = out_temp.I_filled_DRSM;
    out(i).I_filled_ADMM = out_temp.I_filled_ADMM;
    out(i).error_DRSM = out_temp.error_DRSM;
    out(i).error_ADMM = out_temp.error_ADMM;
    
    % 显示原始图像
    axes(ha(i));
    imshow(I);  
    
    % 显示采样后的图像
    axes(ha(3 + i));
    imshow(uint8(out(i).I_sampled));  
    
    % 显示 DRSM 填充后的图像
    axes(ha(6 + i));
    imshow(uint8(out(i).I_filled_DRSM));  
    
    % 显示 ADMM 填充后的图像
    axes(ha(9 + i));
    imshow(uint8(out(i).I_filled_ADMM));  
end
