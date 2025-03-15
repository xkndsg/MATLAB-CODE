function ssim_value = compute_ssim(original, recovered)
    if size(original, 3) == 1
        % 灰度图像直接计算 SSIM
        ssim_value = ssim(uint8(recovered), uint8(original));
    else
        % 计算彩色图像的 SSIM，分别计算 R、G、B 频道的 SSIM，取平均值
        ssim_r = ssim(uint8(recovered(:,:,1)), uint8(original(:,:,1)));
        ssim_g = ssim(uint8(recovered(:,:,2)), uint8(original(:,:,2)));
        ssim_b = ssim(uint8(recovered(:,:,3)), uint8(original(:,:,3)));
        ssim_value = mean([ssim_r, ssim_g, ssim_b]);
    end
end
