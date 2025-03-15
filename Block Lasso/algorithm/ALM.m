function [x_sol,out] =ALM(x,A,b,opts)
% min_x  \|x\|_1  s.t.  Ax=b
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
mu = 1/beta;

opts_1.maxit = 28;
opts_1.maxsigma = norm(A*A', 2);

in_iter = 0; res = 1;
tic;
for iter = 1:opts.maxit
    opts_1.eps = 1e-3/(iter^2);
    
    b1 = b + mu*lam;
    [x_new,iter1]= sub_FISTA(x,A,b1,mu,opts_1);
    
    lin = A*x_new  -b;
    lam_new = lam - lin/mu;

    res1 = norm(x-x_new); res2 = norm(lin);
    res = max(res1 , res2);

    x = x_new;  lam = lam_new; 
    out.res(iter) = norm(x-opts.x)/opts.normx;
    out.time(iter) = toc;
    in_iter = in_iter + iter1;
    out.error = res;
    if res <= opts.eps
        break;
    end
    if mod(iter,5)==0 && beta < 1e8
        beta = beta*10;
        mu = mu/10;
    end
end

x_sol =x;
out.iter = iter;
out.iter1 = in_iter;

end

%用于S-ADMM算法中求解子问题
function [x_sol,iter] = sub_FISTA(x,A,b,mu,opts)

t = 1/opts.maxsigma;
gamma = 1;
y = x; 

for iter = 1:opts.maxit
    x_new = prox(y - t * (A'*(A*y - b)) , mu/t);
    gamma_new = (1 + sqrt(1 + 4*gamma^2) ) / 2;
    y = x_new + (gamma-1) / gamma_new * (x_new-x);

    lin = A*x_new-b;
    
    res1 = norm(x-x_new); res2 = norm(lin)/10;
    res = max(res1 , res2);
    
    x = x_new;  gamma = gamma_new;
    
    if res <= opts.eps
        break;
    end
    
end

x_sol = x;

end
