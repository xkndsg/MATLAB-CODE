function out = DRSM(M, index, beta, rho, kn, tau, max_iter, tol)
   % Default parameters
   %  if nargin < 7
   %      tol = 1e-5;
   %  end
   %  if nargin < 6
   %      max_iter = 1000;
   %  end
   %  if nargin < 5
   %      tao = 1.0;
   %  end
   %  if nargin < 4
   %      rho = 1.0;
   %  end
   %  if nargin < 3
   %      beta = 01;
   %  end

    % Initialize variables
    [m, n] = size(M);
    X = M;
    Y = X;
    U = zeros(m, n);
    % Vx = X;
    % Vy = Y;
    Zx = X;
    Zy = U;
    out.iter=1;
    out.time=0;
    sv = 100;
    for k = 1:max_iter
        Zy_prev = Zy;
        Zx_prev = Zx;
        % v-subproblem
        tic;
        Vx = ((beta + rho)/ (2 * beta + rho) ).* X + (beta / (2 * beta + rho)).* Y;
        Vy = (beta / (beta + rho)).* Vx + (rho/ (beta + rho)) .* Y;
        Vx = Vx.*~index+M.*index;
        % Vy = Vy.*~index+M.*index;
        % z-subproblem
        Z1 = 2 * Vy - Y;
        Z2 = 2 * Vx - X;
        % Z2 = Z2.*~index+M.*index;
        out1 = soft_thresholding(Z2, tau/rho,sv);
        Zx = out1.X;
        sv = out1.svp;
        t=toc;
        out.time = out.time+t;
        Zy = project_to_observed(Z1, M, index);          
        % Update u
        X = X + kn*(Zx - Vx);
        tic;
        Y = Y + kn*(Zy - Vy);
        % X = X.*~index+M.*index;
        % Y = Y.*~index+M.*index;

        e1 = norm((Zx -Zx_prev),"fro") / (norm(Zx,"fro"));
        e2 = norm((Zy -Zy_prev),"fro")/(norm(Zy,"fro"));
        % e1 = norm(~index.*(Zx -Zx_prev),"fro") / (norm(~index.*Zx,"fro"));
        % e2 = norm(~index.*(Zy -Zy_prev),"fro")/(norm(~index.*Zy,"fro"));
        % e1 = norm(~index.*(Zx -Zx_prev),"fro") / (1+norm(~index.*Zx,"fro"));
        % e2 = norm(~index.*(Zy -Zy_prev),"fro")/(1+norm(~index.*Zy,"fro"));
        % e1 = norm(index.*(Zx -Zx_prev),"fro") / max(norm(index.*(Zx),"fro"),10);
        % e2 = norm(index.*(Zy -Zy_prev),"fro")/(10*norm([index.*Zx,index.*Zy],"fro"));
        out.err(k) = max(e1,e2);
        t=toc;
        out.time = out.time+t;
        % Check for convergence
        if out.err(k) <= tol || k>=max_iter
            out.iter=k;
            % out.Vx = Vx;
            % out.Vy = Vy;
            out.Zx = Zx;
            out.Zy = Zy;
            out.X = X;
            out.Y = Y;
            % out.iter = k;
            break;
        end
    end
end

