function [xupd, Pupd, RMS, RMS1] = KFTCM_without_classification(A,Bd,C,Q_true,R_true,Xtrue,Ytrue,mc, RMS, RMS1)
% KFTCM
% 使用真实方差矩阵
xupd(:,1) = zeros(4,1);
Pupd(:,:,1) = zeros(4);%eye(4);%zeros(4); %对称就行

%% Stretch 1, A
for k = 2:1501 % A
    xpre(:,k) = A * xupd(:,k-1);
    Ppre(:,:,k) = A * Pupd(:,:,k-1) *A' + Bd*Q_true(:,k-1)*Bd'; %是没有路面识别，只知道第一段的标称值，不知道后面两段的标称值

    S(:,:,k) = C*Ppre(:,:,k)*C' + diag(R_true(:,k));
    K(:,:,k) = Ppre(:,:,k)*C' / S(:,:,k);

    xupd(:,k) = xpre(:,k) + K(:,:,k)*( Ytrue(:, k) - C*xpre(:,k) );
    Pupd(:,:,k) = Ppre(:,:,k) - K(:,:,k)*C*Ppre(:,:,k); %

    RMS(1, k, mc) = sum((xupd(1:2, k) -Xtrue(1:2, k)).^2,1);
    RMS1(1, k, mc) = sum((xupd(3:4, k) -Xtrue(3:4, k)).^2,1);
end

% %% Stretch 2, C
% for k = 502:1001 % A
%     xpre(:,k) = A * xupd(:,k-1);
%     Ppre(:,:,k) = A * Pupd(:,:,k-1) *A' + Bd*Q_true(:,k-1)*Bd';
% 
%     S(:,:,k) = C*Ppre(:,:,k)*C' + diag(R_true(:,k));
%     K(:,:,k) = Ppre(:,:,k)*C' / S(:,:,k);
% 
%     xupd(:,k) = xpre(:,k) + K(:,:,k)*( Ytrue(:, k) - C*xpre(:,k) );
%     Pupd(:,:,k) = Ppre(:,:,k) - K(:,:,k)*C*Ppre(:,:,k); %
% 
%     RMS(1, k, mc) = sum((xupd(1:2, k) -Xtrue(1:2, k)).^2, 1);
%     RMS1(1, k, mc) = sum((xupd(3:4, k) -Xtrue(3:4, k)).^2, 1);
% end
% 
% %% Stretch 3, E
% for k = 1002:1501 % A
%     xpre(:,k) = A * xupd(:,k-1);
%     Ppre(:,:,k) = A * Pupd(:,:,k-1) *A' + Bd*Q_true(:,k-1)*Bd';
% 
%     S(:,:,k) = C*Ppre(:,:,k)*C' + diag(R_true(:,k));
%     K(:,:,k) = Ppre(:,:,k)*C' / S(:,:,k);
% 
%     xupd(:,k) = xpre(:,k) + K(:,:,k)*( Ytrue(:, k) - C*xpre(:,k) );
%     Pupd(:,:,k) = Ppre(:,:,k) - K(:,:,k)*C*Ppre(:,:,k); %
% 
%     RMS(1, k, mc) = sum((xupd(1:2, k) -Xtrue(1:2, k)).^2, 1);
%     RMS1(1, k, mc) = sum((xupd(3:4, k) -Xtrue(3:4, k)).^2, 1);
% end