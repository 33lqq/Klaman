function  [xupd_3, Pupd_3, R_hat3, RMS_3, RMS_3_1, SRNFN_P3, SRNFN_R3] = VBAKF_R_without_classification(A,Bd,C,Q,R,Q_cons,R_cons,R_true,Xtrue,Ytrue,Pupd,mc,  RMS_3, RMS_3_1, SRNFN_P3, SRNFN_R3)
% parameter;
nx = size(Xtrue,1);
ny = size(Ytrue,1);

zou = 1-exp(-4);%0.9;
tao = 0.8; %0.8
nami = 0.98;

Nsim = 6;
%% ini, A ,既然没有识别了，也就不需要分段了
xupd_3(:,:,1) = zeros(4,1);
Pupd_3(:,:,1) = zeros(4); %怎么初始化？

uupd_3(1) = 5; %3
Uupd_3(:,:,1) = diag(R(:,5));

R_hat_3(:,:,1) = diag(R(:,5));
for k = 2:1501 % epoch 500
    xpre_3(:,:,k) = A*xupd_3(:,:,k-1);
    Ppre_3(:,:,k) = A*Pupd_3(:,:,k-1)*A' + Bd*Q_cons(:,k-1)*Bd';

    upre_3(k) = zou*uupd_3(k-1) ;
    Upre_3(:,:,k)= sqrt(zou*eye(3))*Uupd_3(:,:,k-1)*sqrt(zou*eye(3))';

    uupd_i_3(1) = upre_3(k);
    Uupd_i_3(:,:,1) = Upre_3(:,:,k);
    R_i_3(:,:,1) = (uupd_i_3(1)-ny-1)\Uupd_i_3(:,:,1);

    for i=1:Nsim  % Ppre_i(:,:,i),Ppre(:,:,k)
        S_3(:,:,i) = C*Ppre_3(:,:,k)*C' + R_i_3(:,:,i); % 主要是这个引起的发散
        K_3(:,:,i) = Ppre_3(:,:,k)*C' / ( C*Ppre_3(:,:,k)*C' + R_i_3(:,:,i) );

        xupd_i_3(:,:,i) = xpre_3(:,:,k) + K_3(:,:,i)*(Ytrue(:, k) - C*xpre_3(:,:,k) );
        Pupd_i_3(:,:,i) = Ppre_3(:,:,k) - K_3(:,:,i)*C*Ppre_3(:,:,k); %

        % update: R_i
        B_i_3(:,:,i+1) = C*Pupd_i_3(:,:,i)*C' + (Ytrue(:, k) - C*xupd_i_3(:,:,i))*...
            (Ytrue(:, k) - C*xupd_i_3(:,:,i))';
        uupd_i_3(i+1) = upre_3(k) + 1;
        Uupd_i_3(:,:,i+1) = B_i_3(:,:,i+1) + Upre_3(:,:,k);
        R_i_3(:,:,i+1) = Uupd_i_3(:,:,i+1)/(uupd_i_3(i+1) - ny - 1);
    end
    xupd_3(:,:,k) = xupd_i_3(:,:,Nsim);
    Pupd_3(:,:,k) = Pupd_i_3(:,:,Nsim);
    uupd_3(k) = uupd_i_3(Nsim+1);Uupd_3(:,:,k) = Uupd_i_3(:,:,Nsim+1);
    R_hat_3(:,:,k) = (uupd_3(k)-ny-1)\Uupd_3(:,:,k); R_hat3(:,k) = diag(R_hat_3(:,:,k),0);
    R_hat3(1,k) =  R_hat3(2,k);

    RMS_3(1, k, mc) = sum((xupd_3(1:2, k) -Xtrue(1:2, k)).^2,1);
    RMS_3_1(1, k, mc) = sum((xupd_3(3:4, k) -Xtrue(3:4, k)).^2,1);

    SRNFN_P3(1, k, mc) = 1/4^2*trace(( Pupd_3(:,:,k) - Pupd(:,:,k) )*( Pupd_3(:,:,k)  - Pupd(:,:,k) )');
    SRNFN_R3(1, k, mc) = 1/3^2*trace(( diag(R_hat3(:,k)) - diag(R_true(:, k)))*( diag(R_hat3(:,k)) - diag(R_true(:, k)))');

end

% %% ini, C , L+1
% uupd_3(501) = 5; %3
% Uupd_3(:,:,501) = diag(R(:,5));
% 
% R_hat_3(:,:,501) = diag(R(:,5));
% for k = 502:1001 %epoch k
%     xpre_3(:,:,k) = A * xupd_3(:,:,k-1);
%     Ppre_3(:,:,k) = A * Pupd_3(:,:,k-1) *A' + Bd*Q_cons(:,k-1)*Bd';
% 
%     upre_3(k) = zou*uupd_3(k-1) ;
%     Upre_3(:,:,k)= sqrt(zou*eye(3))*Uupd_3(:,:,k-1)*sqrt(zou*eye(3))';
% 
%     uupd_i_3(1) = upre_3(k);
%     Uupd_i_3(:,:,1) = Upre_3(:,:,k);
%     R_i_3(:,:,1) = (uupd_i_3(1)-ny-1)\Uupd_i_3(:,:,1);
% 
%     for i=1:Nsim  % Ppre_i(:,:,i),Ppre(:,:,k)
%         S_3(:,:,i) = C*Ppre_3(:,:,k)*C' + R_i_3(:,:,i); % 主要是这个引起的发散
%         K_3(:,:,i) = Ppre_3(:,:,k)*C' / ( C*Ppre_3(:,:,k)*C' + R_i_3(:,:,i) );
% 
%         xupd_i_3(:,:,i) = xpre_3(:,:,k) + K_3(:,:,i)*(Ytrue(:, k) - C*xpre_3(:,:,k) );
%         Pupd_i_3(:,:,i) = Ppre_3(:,:,k) - K_3(:,:,i)*C*Ppre_3(:,:,k); %
% 
%         % update: R_i
%         B_i_3(:,:,i+1) = C*Pupd_i_3(:,:,i)*C' + (Ytrue(:, k) - C*xupd_i_3(:,:,i))*...
%             (Ytrue(:, k) - C*xupd_i_3(:,:,i))';
%         uupd_i_3(i+1) = upre_3(k) + 1;
%         Uupd_i_3(:,:,i+1) = B_i_3(:,:,i+1) + Upre_3(:,:,k);
%         R_i_3(:,:,i+1) = Uupd_i_3(:,:,i+1)/(uupd_i_3(i+1) - ny - 1);
%     end
%     xupd_3(:,:,k) = xupd_i_3(:,:,Nsim);
%     Pupd_3(:,:,k) = Pupd_i_3(:,:,Nsim);
%     uupd_3(k) = uupd_i_3(Nsim+1);Uupd_3(:,:,k) = Uupd_i_3(:,:,Nsim+1);
%     R_hat_3(:,:,k) = (uupd_3(k)-ny-1)\Uupd_3(:,:,k);R_hat3(:,k) = diag(R_hat_3(:,:,k),0);
%     R_hat3(1,k) =  R_hat3(2,k);
% 
%     RMS_3(1, k, mc) = sum((xupd_3(1:2, k) -Xtrue(1:2, k)).^2,1);
%     RMS_3_1(1, k, mc) = sum((xupd_3(3:4, k) -Xtrue(3:4, k)).^2,1);
% 
%     SRNFN_P3(1, k, mc) = 1/4^2*trace(( Pupd_3(:,:,k)  - Pupd(:,:,k) )*( Pupd_3(:,:,k)  - Pupd(:,:,k) )');
%     SRNFN_R3(1, k, mc) = 1/3^2*trace(( diag(R_hat3(:,k)) - diag(R_true(:, k)))*( diag(R_hat3(:,k)) - diag(R_true(:, k)))');
% 
% end
% %% ini, E , 2*L+1
% uupd_3(1001) = 5; %3
% Uupd_3(:,:,1001) = diag(R(:,5));
% 
% R_hat_3(:,:,1001) = diag(R(:,5));
% for k = 1002:1501 %epoch k
%     xpre_3(:,:,k) = A * xupd_3(:,:,k-1);
%     Ppre_3(:,:,k) = A * Pupd_3(:,:,k-1) *A' + Bd*Q_cons(:,k-1)*Bd';
%     % variational measurement update
%     upre_3(k) = zou*uupd_3(k-1) ;
%     Upre_3(:,:,k)= sqrt(zou*eye(3))*Uupd_3(:,:,k-1)*sqrt(zou*eye(3))';
% 
%     uupd_i_3(1) = upre_3(k);
%     Uupd_i_3(:,:,1) = Upre_3(:,:,k);
%     R_i_3(:,:,1) = (uupd_i_3(1)-ny-1)\Uupd_i_3(:,:,1);
% 
%     for i=1:Nsim  % Ppre_i(:,:,i),Ppre(:,:,k)
%         S_3(:,:,i) = C*Ppre_3(:,:,k)*C' + R_i_3(:,:,i); % 主要是这个引起的发散
%         K_3(:,:,i) = Ppre_3(:,:,k)*C' / ( C*Ppre_3(:,:,k)*C' + R_i_3(:,:,i) );
% 
%         xupd_i_3(:,:,i) = xpre_3(:,:,k) + K_3(:,:,i)*(Ytrue(:, k) - C*xpre_3(:,:,k) );
%         Pupd_i_3(:,:,i) = Ppre_3(:,:,k) - K_3(:,:,i)*C*Ppre_3(:,:,k); %
% 
%         % update: R_i
%         B_i_3(:,:,i+1) = C*Pupd_i_3(:,:,i)*C' + (Ytrue(:, k) - C*xupd_i_3(:,:,i))*...
%             (Ytrue(:, k) - C*xupd_i_3(:,:,i))';
%         uupd_i_3(i+1) = upre_3(k) + 1;
%         Uupd_i_3(:,:,i+1) = B_i_3(:,:,i+1) + Upre_3(:,:,k);
%         R_i_3(:,:,i+1) = Uupd_i_3(:,:,i+1)/(uupd_i_3(i+1) - ny - 1);
%     end
%     xupd_3(:,:,k) = xupd_i_3(:,:,Nsim);
%     Pupd_3(:,:,k) = Pupd_i_3(:,:,Nsim);
%     uupd_3(k) = uupd_i_3(Nsim+1);Uupd_3(:,:,k) = Uupd_i_3(:,:,Nsim+1);
%     R_hat_3(:,:,k) = (uupd_3(k)-ny-1)\Uupd_3(:,:,k);  R_hat3(:,k) = diag(R_hat_3(:,:,k),0);
%     R_hat3(1,k) =  R_hat3(2,k);
% 
%     RMS_3(1, k, mc) = sum((xupd_3(1:2, k) -Xtrue(1:2, k)).^2,1);
%     RMS_3_1(1, k, mc) = sum((xupd_3(3:4, k) -Xtrue(3:4, k)).^2,1);
% 
%     SRNFN_P3(1, k, mc) = 1/4^2*trace(( Pupd_3(:,:,k)  - Pupd(:,:,k) )*( Pupd_3(:,:,k)  - Pupd(:,:,k) )');
%     SRNFN_R3(1, k, mc) = 1/3^2*trace(( diag(R_hat3(:,k)) - diag(R_true(:, k)))*( diag(R_hat3(:,k)) - diag(R_true(:, k)))');
% 
% end
