function Y = project_to_observed(Y, M, index)
    % 投影到观测值上，保持观测数据不变
       Y = Y.*~index+M.*index;
end