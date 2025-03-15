%多种算法求解LASSO问题对比
%% parameter settings
clear;clc;
addpath ./algorithm
addpath ./function

for i = 2:2
    for  t = 5:5
        M = 512*t; N = 2048*t;    % matrix dimension M-by-N
        S= 10/64*M;                 % sparsity稀疏度

        for j = 1:10
            %% construct sensing matrix
            A = randn(M,N); A = A / norm(A); % Gaussian matrix 
            x_true= zeros(N,1); xs = randn(S,1); p = randi(N,S,1);
            x_true(p) = xs;  %true vector
            b = A*x_true;
            %kk = 1/norm(b,'inf');x_true = x_true*kk; b = b*kk;
            switch i
                case 1
                    eps = 1e-12; x_star = randn(N, 1);
                case 2
                    eps = 1e-8; x_star = zeros(N, 1);
            end
           %% 参数
            opts.x = x_true;
            opts.normx = norm(x_true);
            opts.maxit = 1000; 
            opts.eps = eps; 
            opts.tol = eps; 

            %%% ALM **********
            [~, out_1] = ALM(x_star,A, b, opts);
            [~, out_2] = FISTA(x_star,A, b, opts);

            %%% ADMM **********
            opts.x0 = x_star;
%           [x_3, out_3] = P_ADMM(x_star,A, b, opts);
            opts.rho = 1e-3;
            [~, out_4] = yall1(A, b, opts);
            %%% R-ADMM *********
            [~,out_5] = R_ADMM(x_star,A,b,opts);
            %%% DRSM *********
            opts.kn = 1.618;
            opts.beta =1000;
            opts.rho =1.618;
            [~,out_6] =DRSM(x_star,A,b,opts);
            clear opts
            % print
            zz1(j,:) = [out_1.iter,out_2.iter,out_4.iter,out_5.iter,out_6.iter];%
            zz2(j,:) = [out_1.time(end),out_2.time(end),out_4.time(end),out_5.time(end),out_6.time(end)];%
            disp([t,j]);
        end
        %% p2 rader plot
        load matlab.mat
        zz1 = zz1./2.5;
        z2=zz2.*20;
        X=floor([max(zz1);mean(zz1);min(zz1);max(z2);mean(z2);min(z2)]');
        RC=radarChart(X);
        RC.RLim=[0,350];
        RC.RTick=[0:100:300];
        RC.PropName={'iter-max','iter-aver','iter-min','time-max','time-aver','time-min'};
        RC.ClassName={'ALM','FISTA','ADMM','R-ADMM','DRSM'};
        RC=RC.draw();
        RC.legend();
        colorList=['gkcrb'];
        markerList=['x*^odv'];
        for n=1:RC.ClassNum
            RC.setPatchN(n,'Color',colorList(n),'MarkerFaceColor',colorList(n));
            RC.setPatchN(n,'Marker',markerList(n));
        end
% %
%         数据可视化
%         STEPSIZE = 10;
%         figure;
%         set(gcf,'position',[403 119 640 512])
%         X_range = get_range(length(out_1.res),STEPSIZE);
%         semilogy(X_range,out_1.res(X_range),'-yx', 'LineWidth',1.5,'MarkerSize',8);
%         hold on
%         X_range = get_range(length(out_2.res),STEPSIZE);
%         semilogy(X_range,out_2.res(X_range),'-k*', 'LineWidth',1.5,'MarkerSize',8);
%         hold on
%         X_range = get_range(length(out_3.res),STEPSIZE);
%         semilogy(X_range,out_3.res(X_range),'-g^', 'LineWidth',1.5,'MarkerSize',8);
%         hold on
%         X_range = get_range(length(out_4.res),STEPSIZE);
%         semilogy(X_range,out_4.res(X_range),'-co', 'LineWidth',1.5,'MarkerSize',8);
%         hold on
%         X_range = get_range(length(out_5.res),STEPSIZE);
%         semilogy(X_range,out_5.res(X_range),'-rd', 'LineWidth',1.5,'MarkerSize',8);
%         hold on
%         X_range = get_range(length(out_6.res),STEPSIZE);
%         semilogy(X_range,out_6.res(X_range),'-bv', 'LineWidth',1.5,'MarkerSize',8);
%         hold on
%         xlabel('iterations','fontsize',16)
%         ylabel('Primal relative error $\|x^k-x^*\|/\|x^*\|$.','interpreter','latex','fontsize',16);
%         h = legend('ALM','FISTA','ADMM','R-ADMM','DRSM','location','northeast','interpreter','latex','fontsize',12);
%         str1 = strcat('(m,n,s)=(',num2str(M),'*',num2str(N),',',num2str(S),')');
%         title(str1);
%%
        %p1 error plot
        load matlab.mat
        STEPSIZE = 10;
        figure;
        set(gcf,'position',[403 119 640 512])
        X_range = get_range(length(out_1.res),STEPSIZE);
        semilogy(out_1.time(X_range),out_1.res(X_range),'-g^', 'LineWidth',2,'MarkerSize',8);
        hold on
        X_range = get_range(length(out_2.res),STEPSIZE);
        semilogy(out_2.time(X_range),out_2.res(X_range),'-kx', 'LineWidth',2,'MarkerSize',8);
        hold on
%         X_range = get_range(length(out_3.res),STEPSIZE);
%         semilogy(out_3.time(X_range),out_3.res(X_range),'-g^', 'LineWidth',1.5,'MarkerSize',8);
%         hold on
        X_range = get_range(length(out_4.res),STEPSIZE);
        semilogy(out_4.time(X_range),out_4.res(X_range),'-co', 'LineWidth',2,'MarkerSize',8);
        hold on
        X_range = get_range(length(out_5.res),STEPSIZE);
        semilogy(out_5.time(X_range),out_5.res(X_range),'-rd', 'LineWidth',2,'MarkerSize',8);
        hold on
        X_range = get_range(length(out_6.res),STEPSIZE);
        semilogy(out_6.time(X_range),out_6.res(X_range),'-bv', 'LineWidth',2,'MarkerSize',8);
        hold on
        xlabel('Computing Time (s)','fontsize',16)
        ylabel('Primal relative error $\|x^k-x^*\|/\|x^*\|$.','interpreter','latex','fontsize',16);
        h = legend('ALM','FISTA','ADMM','R-ADMM','DRSM','location','northeast','interpreter','latex','fontsize',12);
        str1 = strcat('(m,n,s)=(',num2str(M),',',num2str(N),',',num2str(S),')');
        title(str1);
    end
end

