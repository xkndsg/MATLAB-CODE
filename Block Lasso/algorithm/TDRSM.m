function [x_sol, out] = TDRSM(x, A, b, opts)
    if ~isfield(opts, 'maxit'); opts.maxit = 1000; end
    if ~isfield(opts, 'eps'); opts.eps = 1e-4; end
    if ~isfield(opts, 'ratio'); opts.ratio = 1/2; end
    if ~isfield(opts, 'beta'); opts.beta = 1; end
    if ~isfield(opts, 'rho'); opts.rho = 0.5; end
    if ~isfield(opts, 'kn'); opts.kn = 0.5; end

    [~, n] = size(A);
    n1 = floor(opts.ratio * n);
    xn = x(1:n1); yn = x(n1+1:end);
    zx = xn; zy = yn;
    AA = A(:, 1:n1); AB = A(:, n1+1:end);
    beta = opts.beta;
    rho = opts.rho;
    kn = opts.kn;

    S = @(x, tau) sign(x) .* max(abs(x) - tau, 0);

    inv_term = inv(beta * (AB' * AB) + rho * eye(size(AB, 2)));
    N_inv = inv(beta * (AA' * AA) + rho * eye(size(AA, 2)) - beta^2 * AA' * AB * inv_term * AB' * AA);

    tic;
    for iter = 1:opts.maxit
        zx1 = zx; zy1 = zy;
        % Solve for z-subproblem with proximal operator
        zx = S(xn, 1 / rho);
        zy = S(yn, 1 / rho);
        
        % Update xtilda
        N1 = beta * (AA' * b) + rho * (2*zx - xn);
        N2 = beta * (AB' * b) + rho *(2*zy- yn);
        xtilda = N1 - beta * AA' * (AB * (inv_term * N2));
        
        % Solve for v-subproblem
        vx = N_inv * xtilda;
        vy = -beta * (inv_term * (AB' * (AA * vx))) + inv_term * N2;
        
        % u-subproblem update
        xn = xn + kn * (vx - zx);
        yn = yn + kn * (vy - zy);

        % Calculate error
        x_sol = [zx; zy];
        e1 = norm(x_sol - [zx1; zy1]) / max(norm(x_sol), 1);
        e2 = norm(A * x_sol - b);
        e = max(e1, e2);
        out.error(iter )=e;
        out.time(iter) = toc;
        out.res(iter) = norm(x_sol - opts.x) / opts.normx;

        if e <= opts.eps
            break;
        end
    end
    out.iter = iter;
end
