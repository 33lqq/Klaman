function [xupd_1, Pupd_1, RMS_1, RMS_1_1, SRNFN_P1, SRNFN_R1] = KFNCM(A,Bd,C,Q_cons,R_cons,R_true,Xtrue,Ytrue,Pupd,N,mc, RMS_1, RMS_1_1, SRNFN_P1, SRNFN_R1)
%KFNCM
%使用固定的方差矩阵
xupd_1(:,1) = zeros(4,1);
Pupd_1(:,:,1) = zeros(4);%eye(4);%zeros(4); %对称就行
%% Stretch 1, A
for k = 2:N % A
    xpre_1(:,k) = A * xupd_1(:,k-1);
    Ppre_1(:,:,k) = A * Pupd_1(:,:,k-1) *A' + Bd*Q_cons(:,k-1)*Bd';

    S_1(:,:,k)  = C*Ppre_1(:,:,k)*C' + diag(R_cons(:,k));
    K_1(:,:,k) = Ppre_1(:,:,k)*C' / S_1(:,:,k);

    xupd_1(:,k) = xpre_1(:,k) + K_1(:,:,k)*(Ytrue(:, k) - C*xpre_1(:,k) );
    Pupd_1(:,:,k) = Ppre_1(:,:,k) - K_1(:,:,k)*C*Ppre_1(:,:,k); %

    RMS_1(1, k, mc) = sum((xupd_1(1:2, k) -Xtrue(1:2, k)).^2,1);
    RMS_1_1(1, k, mc) = sum((xupd_1(3:4, k) -Xtrue(3:4, k)).^2,1);

    SRNFN_P1(1, k, mc) = 1/4^2*trace( ( Pupd_1(:,:,k) - Pupd(:,:,k) )*( Pupd_1(:,:,k)  - Pupd(:,:,k) )');
    SRNFN_R1(1, k, mc) = 1/3^2*trace( ( diag(R_cons(:,k)) - diag(R_true(:, k)))*( diag(R_cons(:,k)) - diag(R_true(:, k)))');
end

%% Stretch 2, C
% for k = 502:1001 % A
%     xpre_1(:,k) = A * xupd_1(:,k-1);
%     Ppre_1(:,:,k) = A * Pupd_1(:,:,k-1) *A' + Bd*Q_cons(:,k-1)*Bd';
% 
%     S_1(:,:,k)  = C*Ppre_1(:,:,k)*C' + diag(R_cons(:,k));
%     K_1(:,:,k) = Ppre_1(:,:,k)*C' / S_1(:,:,k);
% 
%     xupd_1(:,k) = xpre_1(:,k) + K_1(:,:,k)*(Ytrue(:, k) - C*xpre_1(:,k) );
%     Pupd_1(:,:,k) = Ppre_1(:,:,k) - K_1(:,:,k)*C*Ppre_1(:,:,k); %
% 
%     RMS_1(1, k, mc) = sum((xupd_1(1:2, k) -Xtrue(1:2, k)).^2, 1);
%     RMS_1_1(1, k, mc) = sum((xupd_1(3:4, k) -Xtrue(3:4, k)).^2, 1);
% 
%     SRNFN_P1(1, k, mc) = 1/4^2*trace(( Pupd_1(:,:,k)  - Pupd(:,:,k) )*( Pupd_1(:,:,k)  - Pupd(:,:,k) )');
%     SRNFN_R1(1, k, mc) = 1/3^2*trace(( diag(R_cons(:,k)) - diag(R_true(:, k)))*( diag(R_cons(:,k)) - diag(R_true(:, k)))');
% end
% 
% %% Stretch 3, E
% for k = 1002:1501 % A
%     xpre_1(:,k) = A * xupd_1(:,k-1);
%     Ppre_1(:,:,k) = A * Pupd_1(:,:,k-1) *A' + Bd*Q_cons(:,k-1)*Bd';
% 
%     S_1(:,:,k)  = C*Ppre_1(:,:,k)*C' + diag(R_cons(:,k));
%     K_1(:,:,k) = Ppre_1(:,:,k)*C' / S_1(:,:,k);
% 
%     xupd_1(:,k) = xpre_1(:,k) + K_1(:,:,k)*(Ytrue(:, k) - C*xpre_1(:,k) );
%     Pupd_1(:,:,k) = Ppre_1(:,:,k) - K_1(:,:,k)*C*Ppre_1(:,:,k); %
% 
%     RMS_1(1, k, mc) = sum((xupd_1(1:2, k) -Xtrue(1:2, k)).^2, 1);
%     RMS_1_1(1, k, mc) = sum((xupd_1(3:4, k) -Xtrue(3:4, k)).^2, 1);
% 
%     SRNFN_P1(1, k, mc) = 1/4^2*trace(( Pupd_1(:,:,k)  - Pupd(:,:,k) )*( Pupd_1(:,:,k)  - Pupd(:,:,k) )');
%     SRNFN_R1(1, k, mc) = 1/3^2*trace(( diag(R_cons(:,k)) - diag(R_true(:, k)))*( diag(R_cons(:,k)) - diag(R_true(:, k)))');
% end