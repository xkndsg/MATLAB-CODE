function [x_sol, out] = DRSM(x, A, b, opts)
    % if ~isfield(opts, 'maxit'); opts.maxit = 1000; end
    % if ~isfield(opts, 'eps'); opts.eps = 1e-4; end
    if ~isfield(opts, 'ratio'); opts.ratio = 1/2; end
    % if ~isfield(opts, 'beta'); opts.beta = 100; end
    % if ~isfield(opts, 'rho'); opts.rho = 0.8; end
    % if ~isfield(opts, 'kn'); opts.kn = 1; end
    [~, n] = size(A);
    n1 = floor(opts.ratio * n);
    xn = x(1:n1); yn = x(n1+1:end);
    zx = xn; zy = yn;
    A1 = A(:, 1:n1); A2 = A(:, n1+1:end);
    m1 = size(A1,2);m2 = size(A2,2);
    beta = opts.beta;
    rho = opts.rho; 
    br = beta/rho;
    kn = opts.kn;
    A1tA1 = A1'*A1;
    % A2A2t = A2*A2';
    A2tA2 = A2'*A2;
    A1tA2 = A1'*A2;
    A2tA1 = A2'*A1;
    A2tb = A2'*b;
    A1tb = A1'*b;
    brA2tA1 =-br.*A2tA1; 
    brA1tb = br .* A1tb;
    brA2tb = br .* A2tb;
    br2A2tb = br^2 .* A2tb;
    S = @(x, tau) sign(x) .* max(abs(x) - tau, 0);
    % inv_term = inv(br * (A2tA2) + eye(m2));
    inv_term = (br * A2tA2 + eye(m2)) \ eye(m2);
    % N_inv = inv(br*(A1tA1)+eye(m1)-(br^2) * A1tA2 *inv_term *A2tA1);
    N_inv = (br*(A1tA1)+eye(m1)-(br^2) * A1tA2 *inv_term *A2tA1)\eye(m1);
    M1 = A1tA2*inv_term * br2A2tb;
    M2 = A1tA2 * inv_term ;
    M3 = inv_term * brA2tA1;
    M4 = inv_term * brA2tb;
    time = 0; 
    for iter = 1:opts.maxit
        zx1 = zx; zy1 = zy;
        tic;
        xtilda =  brA1tb +  xn - M1- M2* ( br .* yn);
        % Solve for v-subproblem
        vx = N_inv * xtilda;
        vy = M3* vx + M4 + inv_term * yn;
        % Solve for z-subproblem 
        zx = S( 2 .* vx - xn, 1 / rho);
        time = time +toc;
        zy = S(2 .* vy - yn, 1/ rho);
        % u-subproblem 
        tic;
        xn = xn + kn * (zx - vx);
        time = time +toc;
        yn = yn + kn * (zy - vy);
        % Calculate error
        x_sol = [zx;zy];
        e1 = norm(x_sol - [zx1; zy1]) / max(norm(x_sol), 1);
        tic;
        e2 = norm(A*x_sol-b);
        e = max(e1, e2); 
        out.error(iter) = e;
        out.res(iter) = norm(x_sol - opts.x) / opts.normx;
        time = time +toc;
        out.time(iter) = time;
        if e <= opts.eps
            break;
        end
        if mod(iter,15)==0 && beta < 1e11
            beta = beta*10;
            rho = rho*10;
        end 
    end
    out.iter = iter;
end