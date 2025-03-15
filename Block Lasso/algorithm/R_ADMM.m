function [x_sol,out] =R_ADMM(x,A,b,opts)
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

lam = opts.lam;
beta = opts.beta;
mu = 1/beta;

opts_1.maxit = 21;     opts_2.maxit = 21;
opts_1.maxsigma = norm(A1*A1', 2);   opts_2.maxsigma =norm(A2*A2', 2);

in_iter = 0;res=1;
tic;
for iter = 1:opts.maxit
    opts_1.eps = 1e-3/(iter^2);  opts_2.eps = opts_1.eps;
    
    b1 = b - A2*x2  + mu*lam;
    [x1_new,iter1]= sub_FISTA(x1,A1,b1,mu,opts_1);
    
    b2 = b - A1*x1_new + mu*lam;
    [x2_new,~]= sub_FISTA(x2,A2,b2,mu,opts_2);
    
    lin = A1*x1_new + A2*x2_new -b;
    lam_new = lam - 1.* lin/mu;

    res1 = norm([x1;x2]-[x1_new;x2_new]); res2 = norm(lin);
    res = max(res1 , res2);

    x1 = x1_new; x2 = x2_new; lam = lam_new; 
    out.res(iter) = norm([x1; x2]-opts.x)/opts.normx;
    out.time(iter) = toc;
    in_iter = in_iter + iter1;
    
    if res <= opts.eps
        break;
    end
    if mod(iter,5)==0 && beta < 1e8
        beta = beta*10;
        mu = mu/10;
    end
end
out.error = res;
x_sol =[x1; x2];
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