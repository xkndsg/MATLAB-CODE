%多种算法求解LASSO问题对比
%% parameter settings
clear;
clc;
addpath ./algorithm
addpath ./function
fid = fopen(['res\result','.txt'],'w');
fprintf(fid,'\n');
fprintf(fid,"M = 512*t N = 2048*t K = 10/64*M \n\n");
fprintf(fid,"Method & time & iter & error\n\n");
for t = 1:10
        M = 512*t; N = 2048*t;    % matrix dimension M-by-N
        K = 10/64*M;                 % sparsity
        for j = 1:2
            %%construct sensing matrix
            A = randn(M,N); % Gaussian matrix
            A = A / norm(A);    
            x_true= zeros(N,1); %true vector
            xs = randn(K,1);
            p = randi(N,K,1);
            x_true(p) = xs;
            b = A*x_true;
            kk = 1/norm(b,'inf');
            x_true = x_true*kk; b = b*kk;
            eps = 1e-8; x_star = randn(N, 1);
            
            % parameter
            opts.x = x_true;
            opts.normx = norm(x_true);
            opts.maxit = 1000;
            opts.eps = eps; 
            opts.tol = eps; 
            opts.x0 = x_star;
            %%% FISTA**********
            [~,out] = FISTA(x_star,A, b, opts);
            %%% ALM **********
            [~, out1] = ALM(x_star,A, b, opts);
            % [x_2, out2] = L_ALM(x_star,A, b, opts);
            % fprintf("   ALM & %4.2f & %d & %4.2e & %4.2e\n", out1.time(end),out1.iter,out1.error,out1.res(end));
            
            %%% ADMM **********
            opts.x0 = x_star;
            opts.rho = 1e-3;
            [~, out4] = yall1(A, b, opts);
            % fprintf("%d &ADMM & %4.2f & %d & %4.2e & %4.2e\n",t, out3.time(end),out3.iter,out3.error(end),out3.res(end));
            
            %%% LR-ADMM *********
            [~,out5] = R_ADMM(x_star,A,b,opts);
            % fprintf("   RADMM & %4.2f & %d & %4.2e & %4.2e\n", out5.time(end),out5.iter,out5.error,out5.res(end));
            
            %%% DRSM ***********
            opts.kn = 1.618;
            opts.beta =1000;
            opts.rho =1.618;
            [~,out7] =DRSM(x_star,A,b,opts);
            % fprintf("  DRSM & %4.2f & %d & %4.2e & %4.2e\n\n", out7.time(end),out7.iter,out7.error(end),out7.res(end));
            clear opts
            % print
            zz1(j,:) = [out1.iter,out5.iter,out7.iter,out4.iter,out.iter];
            zz2(j,:) = [out1.time(end),out5.time(end),out7.time(end),out4.time(end),out.time(end)];%
            zz3(j,:) = [out1.res(end),out5.res(end),out7.res(end),out4.res(end),out.res(end)];%
            disp([t,j]);
       end
       z1 = mean(zz1);  z2 = mean(zz2);  z3 = mean(zz3);
       fprintf(fid,"   FISTA & %3.2f & %3.2f  & %4.2e \n", z1(5),z2(5),z3(5));
       fprintf(fid,"     ALM & %3.2f & %3.2f  & %4.2e \n", z1(1),z2(1),z3(1));
       fprintf(fid,"  R_ADMM & %3.2f & %3.2f  & %4.2e \n", z1(2),z2(2),z3(2));
       fprintf(fid,"    DRSM & %3.2f & %3.2f  & %4.2e \n", z1(3),z2(3),z3(3));
       fprintf(fid,"    ADMM & %3.2f & %3.2f  & %4.2e \n\n", z1(4),z2(4),z3(4));
       clear x_7 x_5 x_3 x_1 out7 out5 out3 out1
end

% %% test demo
% clear;
% clc;
% addpath ./algorithm
% addpath ./function
% % randn('state',0);   rand('twister',0);
% t=1;
%         M = 512*t; N = 2048*t;    % matrix dimension M-by-N
%         K = 10/64*M;                 % sparsity
%         % fprintf("M = 512*t N = 2048*t K = 10/64*M\n\n");
%         % fprintf("Method & time & iter & error\n\n");
%             %%construct sensing matrix
%             A = randn(M,N); % Gaussian matrix
%             A = A / norm(A);    
%             x_true= zeros(N,1); %true vector
%             xs = randn(K,1);
%             p = randi(N,K,1);
%             x_true(p) = xs;
%             b = A*x_true;
%             kk = 1/norm(b,'inf');
%             x_true = x_true*kk; b = b*kk;
%             eps = 1e-8; x_star = randn(N, 1);
%             % 
%             opts.x = x_true;
%             opts.normx = norm(x_true);
%             opts.maxit = 1000; 
%             opts.eps = eps; 
%             opts.tol = eps; 
%             opts.x0 = x_star;
%             % [x_5,out5] = P_ADMM(x_star,A,b,opts);
%             % opts.kn = 1.618;
%             opts.beta =100;
%             opts.mu = 0.1;
%             % for r =3
%             % opts.rho =r;
%             [x_7,out7] =FISTA(x_star,A,b,opts);
%             fprintf("   FISTA & %4.2f & %d & %4.2e & %4.2e\n\n", out7.time(end),out7.iter,out7.error(end),out7.res(end));
%             % end
% 
% %          
% figure;
% semilogy(out7.error, 'r',LineWidth=1);
% hold on;
% % semilogy(out5.error, 'r',LineWidth=1);
