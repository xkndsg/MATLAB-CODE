function [x_sol,out] = LR_ADMM(x,A,b,opts)
% min_x  \|x\|_1  s.t.  Ax=b
% is equivalent to
% min_{x_1,x_2}  \|x_1\|_1 + \|x_2\|_1  s.t.  A_1x_1 + A_2x_2 = b
%
% INPUT:
% A --- an m x n matrix
% b --- an m-vector
% x --- an initial iteration point
% opts : a structure with fields:
%      opts.maxit: The maximum number of iterations
%      opts.eps  : tolerance
%      opts.beta : penalty parameter
%      opts.lam  : initial point of Lagrange multiplier
if ~isfield(opts, 'maxit'); opts.maxit = 1000; end
if ~isfield(opts, 'eps'); opts.eps = 1e-4; end
if ~isfield(opts, 'beta'); opts.beta = mean(abs(b)); end
[m,n] = size(A);
if ~isfield(opts, 'lam'); opts.lam = zeros(m,1); end
%
% OUTPUT:
% x_sol: last iterate
% out  : a structure with fields:
%      out.res : The error from the real solution
%      out.time: Time spent
%      out.iter: The number of iteration steps
%
% By Y. Wu

n1 = floor(n/2);
x1 = x(1:n1); x2 = x(n1+1:end);
A1 = A(:, 1:n1);   A2 = A(:, n1+1:end);

gam1 = 1/norm(A1*A1', 2);  gam2 = 1.33/norm(A2*A2', 2);
lam = opts.lam;
beta = opts.beta;
gb1 = gam1/beta; gb2 = gam2/beta;
lin = A1*x1 + A2*x2 -b;

tic;
for iter = 1:opts.maxit

    x1_new = prox(x1 - gam1.*(A1'*(lin - lam/beta)), gb1);
    lin = lin + A1*(x1_new-x1);
    x2_new = prox(x2 - gam2.*(A2'*(lin - lam/beta)), gb2);

    lin = lin + A2*(x2_new-x2);
    lam = lam - beta* lin;

    res1 = norm([x1;x2]-[x1_new;x2_new]); res2 = norm(lin);
    res = max(res1 , res2);
    
    x1 = x1_new; x2 = x2_new;
    out.res(iter) = norm([x1;x2]-opts.x)/opts.normx;
    out.time(iter) = toc;
    
    if res<= opts.eps
        break;
    end
    if mod(iter,26)==0 && beta < 1e8
        gg = 10;
        beta = beta*gg;
        gb1 = gb1/gg; gb2 = gb2/gg;
    end
end

x_sol = [x1;x2];
out.iter = iter;

end


