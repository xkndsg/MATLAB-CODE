function [x_sol, out] = FISTA(x0, A, b,opts)
% Solve LASSO problem:
%   min_x  1/2 * ||Ax - b||_2^2 + mu * ||x||_1
%
% INPUT:
%   x0    - initial guess for x
%   A     - measurement matrix
%   b     - observation vector
%   lambda- regularization parameter
%   opts  - options structure:
%           opts.maxit : maximum iterations
%           opts.eps   : tolerance for convergence
%           opts.tau   : step size (optional, default 1/L)
% OUTPUT:
%   x_sol - the solution
%   out   - output structure:
%           out.res   : residuals
%           out.time  : computation time
%           out.iter  : number of iterations

if ~isfield(opts, 'maxit'); opts.maxit = 1000; end
if ~isfield(opts, 'eps'); opts.eps = 1e-4; end
if ~isfield(opts, 'mu'); opts.mu =   0.1; end
if ~isfield(opts, 'tau') 
% Compute Lipschitz constant L = ||A'*A||_2
 % L = eigs(A'*A, 1); % Lipschitz constant of gradient of f
 L = norm(A'*A,2);
    opts.tau =1 / L; 
end
mu = opts.mu;
tau = opts.tau; % Step size
% [m, n] = size(A);
x_new = x0; % Initialize x
xp = x_new;
% y = x_new; % Initialize y for FISTA
% t = 0.5; % Initialize t for FISTA

tic;
for iter = 1:opts.maxit
    theta = (iter - 1) / (iter + 2);
    y = x_new + theta .* (x_new - xp);
    xp = x_new;
    % % Gradient step
    grad = A' * (A * y - b); 
    x_new = prox(y - tau * grad, tau*mu); % Proximal operator for L1 norm
    % % Update momentum parameter
    % t_new = (1 + sqrt(1 + 4 * t^2)) / 2;
    % % Update y
    % y = x_new + ((t - 1) / t_new) * (x_new - x);
    % Convergence check
    e1 = norm(x_new - xp);
    e2 = norm(A*x_new -b);
    e =max(e1,e2);
    out.res(iter) = norm(x_new-opts.x)/opts.normx;
    out.time(iter) = toc;
    if e < opts.eps
        break;
    end
    % Update variables
    % x = x_new;
    % t = t_new;
    out.error(iter) = e;
     if mod(iter,25)==0 && mu > 1e-11
            mu = mu*0.5;
            % tau = tau*5;
     end 
end

x_sol = x_new; % Final solution
out.iter = iter;

end

function x = prox(z, tau)
% Proximal operator for L1 norm: soft-thresholding
    x = sign(z) .* max(abs(z) - tau, 0);
end
