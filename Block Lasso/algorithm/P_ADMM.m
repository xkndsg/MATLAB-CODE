function [x_sol,out] = P_ADMM(x,A,b,opts)
% min_x  \|x\|_1  s.t.  Au=b  x-u=0
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

u = zeros(n,1);
lam1 = opts.lam;
lam2 = zeros(n,1);

beta1 = 3*opts.beta;
beta2 = opts.beta;

beta = beta1/beta2;
B = (eye(m,m)+beta*A*A');
C = (eye(n,n)-beta*A'/B*A );
tic;

for iter = 1:opts.maxit
    
    x_new = prox(u+lam2/beta2, 1/beta2);
    
    w = beta1*(A'*(b+lam1/beta1)) + beta2*x_new - lam2;
    u_new = C*w./beta2;
    
    lam1_new = lam1 - beta1 * (A*u_new-b);
    lam2_new = lam2 - beta2 * (x_new-u_new);

    res1 = norm(x-x_new)/norm(x);  res2 = norm(A*x_new -b);
    res = max(res1 , res2);
    out.error(iter) = res;
    x = x_new;  u = u_new; 
    lam1 = lam1_new;  lam2 = lam2_new;
    out.res(iter) = norm(x-opts.x)/opts.normx;
    out.time(iter) = toc;
    
    if res<= opts.eps
        break;
    end
    if mod(iter,26)==0 && beta < 10
        gg = 10;
        beta1 = beta1*gg;  beta2 = beta2*gg;
    end
    
end
% out.error = res;
x_sol = x;
out.iter = iter;

end

