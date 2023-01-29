function [xupd_4, Pupd_4, R_hat4, R_hat41, RMS_4, RMS_4_1, RMS_41, RMS_41_1, SRNFN_P4, SRNFN_R4, SRNFN_P41, SRNFN_R41] = My_filter_without_classification(A,Bd,C,Q,R,Q_cons,R_cons,R_true,Xtrue,Ytrue,Pupd,mc, RMS_4, RMS_4_1, RMS_41, RMS_41_1, SRNFN_P4, SRNFN_P41, SRNFN_R4,SRNFN_R41)
% parameter;
nx = size(Xtrue,1);
ny = size(Ytrue,1);

zou = 1-exp(-4);%0.9;
tao = 0.8; %0.8
nami = 0.98;

Nsim = 6;
%% ini, A ，既然没有路面识别了，也就不需要分段了
xupd_4(:,:,1) = zeros(4,1);
Pupd_4(:,:,1) = zeros(4); %怎么初始化？

uupd(1) = 5; %4
Uupd(:,:,1) = diag(R(:,5));


tupd(1) = 6;
Tupd(:,:,1) = A* Pupd_4(:,:,1) *A' + Bd*Q(:,5)*Bd';%mvnrnd(eye(4),eye(4)); %eye(4);


R_hat_4(:,:,1) = diag(R(:,5));
for k = 2:1501 % epoch 500
    xpre_4(:,:,k) = A * xupd_4(:,:,k-1);
%     Ppre_4(:,:,k) = A * Pupd_4(:,:,k-1) *A' + Bd*Q_true(:,k-1)*Bd';
    Ppre_4(:,:,k) = A * Pupd_4(:,:,k-1) *A' + Bd*Q_cons(:,k-1)*Bd';


    upre(k) = zou*uupd(k-1) ;
    Upre(:,:,k)= sqrt(zou*eye(3))*Uupd(:,:,k-1)*sqrt(zou*eye(3))';

    tpre(k) = tao + nx + 1;
    %         tpre(k) = tao*tupd(k-1) ;
    Tpre(:,:,k)= sqrt(tao*eye(4))*Ppre_4(:,:,k)*sqrt(tao*eye(4))';

    uupd_i(1) = upre(k);
    Uupd_i(:,:,1) = Upre(:,:,k);
    R_i(:,:,1) = (uupd_i(1)-ny-1)\Uupd_i(:,:,1);

    tupd_i(1) = tpre(k);
    Tupd_i(:,:,1) = Tpre(:,:,k);
    Ppre_i(:,:,1) =  Ppre_4(:,:,k); %怎么初始？

    for i=1:Nsim  % Ppre_i(:,:,i),Ppre(:,:,k)
        S_4(:,:,i) = C*Ppre_i(:,:,i)*C' + R_i(:,:,i); % 主要是这个引起的发散
        K_4(:,:,i) = Ppre_i(:,:,i)*C' / ( C*Ppre_i(:,:,i)*C' + R_i(:,:,i) );

        xupd_i(:,:,i) = xpre_4(:,:,k) + K_4(:,:,i)*(Ytrue(:, k) - C*xpre_4(:,:,k) );
        Pupd_i(:,:,i) = Ppre_i(:,:,i) - K_4(:,:,i)*C*Ppre_i(:,:,i); %

        % update: P_i
        A_i(:,:,i+1) =  Pupd_i(:,:,i) + (xupd_i(:,:,i) - xpre_4(:,:,k))*...
            (xupd_i(:,:,i) - xpre_4(:,:,k))'; %这个是不对的
        tupd_i(i+1) = tpre(k) + 1;
        Tupd_i(:,:,i+1) = A_i(:,:,i+1) + Tpre(:,:,k);
        Ppre_i(:,:,i+1) = Tupd_i(:,:,i+1)/(tupd_i(i+1) - nx - 1);

        % update: R_i
        B_i(:,:,i+1) = C*Pupd_i(:,:,i)*C' + (Ytrue(:, k) - C*xupd_i(:,:,i))*...
            (Ytrue(:, k) - C*xupd_i(:,:,i))';
        uupd_i(i+1) = upre(k) + 1;
        Uupd_i(:,:,i+1) = B_i(:,:,i+1) + Upre(:,:,k);
        R_i(:,:,i+1) = Uupd_i(:,:,i+1)/(uupd_i(i+1) - ny - 1);
    end
    xupd_4(:,:,k) = xupd_i(:,:,Nsim);
    Pupd_4(:,:,k) = Pupd_i(:,:,Nsim);
    uupd(k) = uupd_i(Nsim+1);Uupd(:,:,k) = Uupd_i(:,:,Nsim+1);
    tupd(k) = tupd_i(Nsim+1);Tupd(:,:,k) = Tupd_i(:,:,Nsim+1);
    R_hat_4(:,:,k) = (uupd(k)-ny-1)\Uupd(:,:,k);
    R_hat4(:,k) = diag(R_hat_4(:,:,k),0);
    R_hat4(1,k) = R_hat4(2,k);

    P_hat_4(:,:,k) = (tupd(k)-nx-1)\Tupd(:,:,k);

    RMS_4(1, k, mc) = sum((xupd_4(1:2, k) -Xtrue(1:2, k)).^2,1);
    RMS_4_1(1, k, mc) = sum((xupd_4(3:4, k) -Xtrue(3:4, k)).^2,1);

    SRNFN_P4(1, k, mc) = 1/4^2*trace(( Pupd_4(:,:,k)  - Pupd(:,:,k) )*( Pupd_4(:,:,k)  - Pupd(:,:,k) )');
    SRNFN_R4(1, k, mc) = 1/3^2*trace(( diag(R_hat4(:,k)) - diag(R_true(:, k)))*( diag(R_hat4(:,k)) - diag(R_true(:, k)))');

end
N = k;
xupd_bwd(:,:,N) =xupd_4(:,:,N); %xpre_k_i{:,Nsmin+1}, xupd_k_i{:,Nsmin+1}
Pupd_bwd(:,:,N) = Pupd_4(:,:,N);%Ppre_k_i{:,Nsmin+1}, Pupd_k_i{:,Nsmin+1}
R_hat_bwd(:,N) = R_hat4(:,N);
for i =N-1:-1:5
    % 矩阵接近奇异值，或者缩放错误。结果可能不准确,矩阵的行列式接近0，接近不可逆
    K_bwd(:,:,i) = Pupd_4(:,:,i)*A'/Ppre_4(:,:,i+1); %R不是后验
    xupd_bwd(:,:,i) = xupd_4(:,:,i) + K_bwd(:,:,i)*(xupd_bwd(:,:,i+1) - A*xupd_4(:,:,i));
    %产生对角线为负的值，原P全是正的
    Pupd_bwd(:,:,i) = Pupd_4(:,:,i) + K_bwd(:,:,i)*(Pupd_bwd(:,:,i+1)...
        - Ppre_4(:,:,i+1))*K_bwd(:,:,i)';

    R_hat_bwd(:,i) = (1-nami)*R_hat4(:,i) + nami*R_hat_bwd(:,i+1);
end
Pupd_4(:,:,5:N) = Pupd_bwd(:,:,5:N);
xupd_4(:,:,5:N) =  xupd_bwd(:,:,5:N);
R_hat41(:,5:N) = R_hat_bwd(:,5:N);
R_hat41(:,1:4) = R_hat4(:,1:4);
% R_hat41(1,1:N)  = R_hat41(2,1:N);
for k=1:N
    RMS_41(1, k, mc) = sum((xupd_4(1:2, k) -Xtrue(1:2, k)).^2,1);
    RMS_41_1(1, k, mc) = sum((xupd_4(3:4, k) -Xtrue(3:4, k)).^2,1);
end
for k=1:N
    SRNFN_P41(1, k, mc) = 1/4^2*trace(( Pupd_4(:,:,k)  - Pupd(:,:,k) )*( Pupd_4(:,:,k)  - Pupd(:,:,k) )');
    SRNFN_R41(1, k, mc) = 1/3^2*trace(( diag(R_hat41(:,k)) - diag(R_true(:, k)))*( diag(R_hat41(:,k)) - diag(R_true(:, k)))');
end

% %% ini, C , L+1
% %     xupd(:,:,L+1) = zeros(4,1);
% %     Pupd(:,:,L+1) = zeros(4); %不用初始化，上一阶段已经计算好了
% 
% uupd(501) = 5; %4
% Uupd(:,:,501) = diag(R(:,5));
% 
% tupd(501) = 6;
% Tupd(:,:,501) = A* Pupd_4(:,:,501) *A' + Bd*diag(Q(:,5))*Bd';%mvnrnd(eye(4),eye(4)); %eye(4);
% 
% R_hat_4(:,:,501) = diag(R(:,5));
% for k = 502:1001 %epoch k
%     xpre_4(:,:,k) = A * xupd_4(:,:,k-1);
%     Ppre_4(:,:,k) = A * Pupd_4(:,:,k-1) *A' + Bd*Q_cons(:,k-1)*Bd';
%     % variational measurement update
%     %         xupd_i(:,:,1) = xpre(:,:,k);
%     %         Pupd_i(:,:,1) = Ppre(:,:,k);
% 
%     upre(k) = zou*uupd(k-1) ;
%     Upre(:,:,k)= sqrt(zou*eye(3))*Uupd(:,:,k-1)*sqrt(zou*eye(3))';
% 
%     tpre(k) = tao + nx + 1;
%     %         tpre(k) = tao*tupd(k-1) ;
%     Tpre(:,:,k)= sqrt(tao*eye(4))*Ppre_4(:,:,k)*sqrt(tao*eye(4))';
% 
%     uupd_i(1) = upre(k);
%     Uupd_i(:,:,1) = Upre(:,:,k);
%     R_i(:,:,1) = (uupd_i(1)-ny-1)\Uupd_i(:,:,1);
% 
%     tupd_i(1) = tpre(k);
%     Tupd_i(:,:,1) = Tpre(:,:,k);
%     Ppre_i(:,:,1) =  Ppre_4(:,:,k); %怎么初始？
% 
%     for i=1:Nsim  % Ppre_i(:,:,i),Ppre(:,:,k)
%         S_4(:,:,i) = C*Ppre_i(:,:,i)*C' + R_i(:,:,i); % 主要是这个引起的发散
%         K_4(:,:,i) = Ppre_i(:,:,i)*C' / ( C*Ppre_i(:,:,i)*C' + R_i(:,:,i) );
% 
%         xupd_i(:,:,i) = xpre_4(:,:,k) + K_4(:,:,i)*(Ytrue(:, k) - C*xpre_4(:,:,k) );
%         Pupd_i(:,:,i) = Ppre_i(:,:,i) - K_4(:,:,i)*C*Ppre_i(:,:,i); %
% 
%         % update: P_i
%         A_i(:,:,i+1) =  Pupd_i(:,:,i) + (xupd_i(:,:,i) - xpre_4(:,:,k))*...
%             (xupd_i(:,:,i) - xpre_4(:,:,k))'; %这个是不对的
%         tupd_i(i+1) = tpre(k) + 1;
%         Tupd_i(:,:,i+1) = A_i(:,:,i+1) + Tpre(:,:,k);
%         Ppre_i(:,:,i+1) = Tupd_i(:,:,i+1)/(tupd_i(i+1) - nx - 1);
% 
%         % update: R_i
%         B_i(:,:,i+1) = C*Pupd_i(:,:,i)*C' + (Ytrue(:, k) - C*xupd_i(:,:,i))*...
%             (Ytrue(:, k) - C*xupd_i(:,:,i))';
%         uupd_i(i+1) = upre(k) + 1;
%         Uupd_i(:,:,i+1) = B_i(:,:,i+1) + Upre(:,:,k);
%         R_i(:,:,i+1) = Uupd_i(:,:,i+1)/(uupd_i(i+1) - ny - 1);
%     end
%     xupd_4(:,:,k) = xupd_i(:,:,Nsim);
%     Pupd_4(:,:,k) = Pupd_i(:,:,Nsim);
%     uupd(k) = uupd_i(Nsim+1);Uupd(:,:,k) = Uupd_i(:,:,Nsim+1);
%     tupd(k) = tupd_i(Nsim+1);Tupd(:,:,k) = Tupd_i(:,:,Nsim+1);
%     R_hat_4(:,:,k) = (uupd(k)-ny-1)\Uupd(:,:,k);R_hat4(:,k) = diag(R_hat_4(:,:,k),0);
%     R_hat4(1,k) = R_hat4(2,k);
%     P_hat_4(:,:,k) = (tupd(k)-nx-1)\Tupd(:,:,k);
% 
%     RMS_4(1, k, mc) = sum((xupd_4(1:2, k) -Xtrue(1:2, k)).^2,1);
%     RMS_4_1(1, k, mc) = sum((xupd_4(3:4, k) -Xtrue(3:4, k)).^2,1);
% 
%     SRNFN_P4(1, k, mc) = 1/4^2*trace(( Pupd_4(:,:,k)  - Pupd(:,:,k) )*( Pupd_4(:,:,k)  - Pupd(:,:,k) )');
%     SRNFN_R4(1, k, mc) = 1/3^2*trace(( diag(R_hat4(:,k)) - diag(R_true(:, k)))*( diag(R_hat4(:,k)) - diag(R_true(:, k)))');
% 
% end
% %% ini, E , 2*L+1
% uupd(1001) = 5; %4
% Uupd(:,:,1001) = diag(R(:,5));
% 
% tupd(1001) = 6;
% Tupd(:,:,1001) = A* Pupd_4(:,:,1001) *A' + Bd*diag(Q(:,1))*Bd';%mvnrnd(eye(4),eye(4)); %eye(4);
% 
% R_hat_4(:,:,1001) = diag(R(:,1));
% for k = 1002:1501 %epoch k
%     xpre_4(:,:,k) = A * xupd_4(:,:,k-1);
%     Ppre_4(:,:,k) = A * Pupd_4(:,:,k-1) *A' + Bd*Q_cons(:,1)*Bd';
%     % variational measurement update
%     upre(k) = zou*uupd(k-1) ;
%     Upre(:,:,k)= sqrt(zou*eye(3))*Uupd(:,:,k-1)*sqrt(zou*eye(3))';
% 
%     tpre(k) = tao + nx + 1;
%     %         tpre(k) = tao*tupd(k-1) ;
%     Tpre(:,:,k)= sqrt(tao*eye(4))*Ppre_4(:,:,k)*sqrt(tao*eye(4))';
% 
%     uupd_i(1) = upre(k);
%     Uupd_i(:,:,1) = Upre(:,:,k);
%     R_i(:,:,1) = (uupd_i(1)-ny-1)\Uupd_i(:,:,1);
% 
%     tupd_i(1) = tpre(k);
%     Tupd_i(:,:,1) = Tpre(:,:,k);
%     Ppre_i(:,:,1) =  Ppre_4(:,:,k); %怎么初始？
% 
%     for i=1:Nsim  % Ppre_i(:,:,i),Ppre(:,:,k)
%         S_4(:,:,i) = C*Ppre_i(:,:,i)*C' + R_i(:,:,i); % 主要是这个引起的发散
%         K_4(:,:,i) = Ppre_i(:,:,i)*C' / ( C*Ppre_i(:,:,i)*C' + R_i(:,:,i) );
% 
%         xupd_i(:,:,i) = xpre_4(:,:,k) + K_4(:,:,i)*(Ytrue(:, k) - C*xpre_4(:,:,k) );
%         Pupd_i(:,:,i) = Ppre_i(:,:,i) - K_4(:,:,i)*C*Ppre_i(:,:,i); %
% 
%         % update: P_i
%         A_i(:,:,i+1) =  Pupd_i(:,:,i) + (xupd_i(:,:,i) - xpre_4(:,:,k))*...
%             (xupd_i(:,:,i) - xpre_4(:,:,k))'; %这个是不对的
%         tupd_i(i+1) = tpre(k) + 1;
%         Tupd_i(:,:,i+1) = A_i(:,:,i+1) + Tpre(:,:,k);
%         Ppre_i(:,:,i+1) = Tupd_i(:,:,i+1)/(tupd_i(i+1) - nx - 1);
% 
%         % update: R_i
%         B_i(:,:,i+1) = C*Pupd_i(:,:,i)*C' + (Ytrue(:, k) - C*xupd_i(:,:,i))*...
%             (Ytrue(:, k) - C*xupd_i(:,:,i))';
%         uupd_i(i+1) = upre(k) + 1;
%         Uupd_i(:,:,i+1) = B_i(:,:,i+1) + Upre(:,:,k);
%         R_i(:,:,i+1) = Uupd_i(:,:,i+1)/(uupd_i(i+1) - ny - 1);
%     end
%     xupd_4(:,:,k) = xupd_i(:,:,Nsim);
%     Pupd_4(:,:,k) = Pupd_i(:,:,Nsim);
%     uupd(k) = uupd_i(Nsim+1);Uupd(:,:,k) = Uupd_i(:,:,Nsim+1);
%     tupd(k) = tupd_i(Nsim+1);Tupd(:,:,k) = Tupd_i(:,:,Nsim+1);
%     R_hat_4(:,:,k) = (uupd(k)-ny-1)\Uupd(:,:,k);  R_hat4(:,k) = diag(R_hat_4(:,:,k),0);
%     R_hat4(1,k) = R_hat4(2,k);
% 
%     P_hat_4(:,:,k) = (tupd(k)-nx-1)\Tupd(:,:,k);
% 
%     RMS_4(1, k, mc) = sum((xupd_4(1:2, k) -Xtrue(1:2, k)).^2,1);
%     RMS_4_1(1, k, mc) = sum((xupd_4(3:4, k) -Xtrue(3:4, k)).^2,1);
% 
%     SRNFN_P4(1, k, mc) = 1/4^2*trace(( Pupd_4(:,:,k)  - Pupd(:,:,k) )*( Pupd_4(:,:,k)  - Pupd(:,:,k) )');
%     SRNFN_R4(1, k, mc) = 1/3^2*trace(( diag(R_hat4(:,k)) - diag(R_true(:, k)))*( diag(R_hat4(:,k)) - diag(R_true(:, k)))');
% 
% end
