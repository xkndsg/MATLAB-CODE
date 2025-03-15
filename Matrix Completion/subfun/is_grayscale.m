function is_gray = is_grayscale(I)
    % 判断是否为灰度图像
    if size(I, 3) == 1
        is_gray = true;
    else
        is_gray = false;
    end
end
