function [x_sol,out] = L_ALM(x,A,b,opts)
% min_x  \|x\|_1  s.t.  Ax=b
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
[m,~] = size(A);
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

lam = opts.lam;
beta = opts.beta;       
gamma = 1/norm(A*A',2);

tic;
for iter = 1:opts.maxit
    
    temp = x - gamma*(A'*(A*x-b-lam/beta));
    x_new = prox(temp,gamma/beta );
    
    lin = A*x_new - b;
    lam= lam - beta* lin;

    res1 = norm(x_new - x)/norm(x);res2 = norm(lin);
    res = max([res1 , res2]);
    
    x = x_new;
    out.res(iter) = norm(x-opts.x)/opts.normx;
    out.time(iter) = toc;
    
    if res <= opts.eps
        break;
    end
    if mod(iter,26)==0 && beta < 1e8
        gg = 10;
        beta = beta*gg;
    end
end

x_sol = x;
out.iter = iter;

end
