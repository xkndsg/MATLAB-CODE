function out = SVT(M,P,T,delta,itermax,tol)
%Single value thresholding algorithm，SVT
% function：solve the following optimization problem
%                  min  ||X||*
%               s.t. P(X-M) = 0
% X: recovered matrix
% M: observed matrix
% T: single value threshold
% delta: step size
% output：X,iterations

% initialization
Y = zeros(size(M));
X = zeros(size(M));
% Y = M;
% X = M;
% 
% if nargin < 3
%     T =  sqrt(n1*n2);
% end
% if nargin < 4
%     delta = 1;
% end
% if nargin < 5
%     itermax = 1000 ;
% end
% if nargin < 6
%     tol = 1e-4;
% end
 % [m, ~] = size(M);
out.time=0;
sv = 100;
for ii = 1:itermax
    X_pre = X;
    Y_pre = Y;
    % singular value threshold operation
    tic;
    % X = soft_thresholding(Y, T,100+r+m/100);
    out1 = soft_thresholding(Y, T,sv);
    X = out1.X;
    sv = out1.svp;
    % opts = struct('MaxIterations', 100);  % 增加最大迭代次数
    % [U, S, V] = svds(Y, 18) ;
    % S = sign(S) .* max(abs(S) - T, 0) ;
    % X = U*S*V' ;
    % update auxiliary matrix Y
    Y = Y + delta* P.* (M-X);
    Y = P.*Y ;
    
    % computer error
    % error= norm( P.* (M-X),'fro' )/norm( P.* M,'fro' );
    e1= norm((X_pre-X),'fro' )/(norm( X,'fro'));
    e2= norm((Y_pre-Y),'fro' )/(norm( Y,'fro'));
    out.X = X;
    out.err(ii) = max(e1,e2);
    t=toc;
    out.time=out.time+t;
    out.iter = ii ;
    if out.err(ii) <tol
      out.iter = ii ;  break;
    end
    % update iterations
    
end
end