%一范数邻近算子
% Shrinkage operator for $h(x)=\mu\|x\|_1$
function y = prox(x, mub)
y = max(abs(x) - mub, 0);
y = sign(x) .* y;
end
