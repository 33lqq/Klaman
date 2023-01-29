function [xupd_2, Pupd_2, R_hat2, RMS_2, RMS_2_1, SRNFN_P2, SRNFN_R2] = RSAKF_without_classification(A,Bd,C,Q_cons,R_cons,R_true,Xtrue,Ytrue,Pupd,mc, RMS_2, RMS_2_1, SRNFN_P2, SRNFN_R2)
m=20;
xupd_2(:,1) = zeros(4,1);
Pupd_2(:,:,1) = zeros(4);%eye(4);%zeros(4); %对称就行
for i=1:m
    R_hat_2(:,:,i) = diag( R_cons(:,5) );
    R_hat2(:,i) = R_cons(:,5);
end
%% Stretch 1, A
for k = 2:1501 % A
    xpre_2(:,k) = A*xupd_2(:,k-1);
    Ppre_2(:,:,k) = A*Pupd_2(:,:,k-1)*A' + Bd*Q_cons(:,k-1)*Bd';

    if k<m
        S_2(:,:,k)  = C*Ppre_2(:,:,k)*C' + diag(R_cons(:,k));
        K_2(:,:,k) = Ppre_2(:,:,k)*C' / S_2(:,:,k);
    else % m
        S_2(:,:,k)  = C*Ppre_2(:,:,k)*C' + R_hat_2(:,:,k);
        K_2(:,:,k) = Ppre_2(:,:,k)*C' / S_2(:,:,k);
    end

    xupd_2(:,k) = xpre_2(:,k) + K_2(:,:,k)*(Ytrue(:, k) - C*xpre_2(:,k) );
    Pupd_2(:,:,k) = Ppre_2(:,:,k) - K_2(:,:,k)*C*Ppre_2(:,:,k); %

    e_y(:, k) =  Ytrue(:, k) -  C*xupd_2(:,k); %k=2才开始有值

    if k>=m
        C_v = 0;
        for j=1:m
            C_v = C_v + e_y(:,k-j+1)*e_y(:,k-j+1)';
        end
        R_hat_2(:,:,k+1) = C_v/m + C*Pupd_2(:,:,k)*C'; % k时刻的估计给k+1时刻使用
        R_hat2(:,k+1) = diag(R_hat_2(:,:,k+1),0);
    end
    RMS_2(1, k, mc) = sum((xupd_2(1:2, k) -Xtrue(1:2, k)).^2,1);
    RMS_2_1(1, k, mc) = sum((xupd_2(3:4, k) -Xtrue(3:4, k)).^2,1);

    SRNFN_P2(1, k, mc) = 1/4^2*trace(( Pupd_2(:,:,k)  - Pupd(:,:,k) )*( Pupd_2(:,:,k)  - Pupd(:,:,k) )');
    SRNFN_R2(1, k, mc) = 1/3^2*trace(( diag(R_hat2(:,k)) - diag(R_true(:, k)))*( diag(R_hat2(:,k)) - diag(R_true(:, k)))');

end

% %% Stretch 2, C
% for k = 502:1001 % A
%     xpre_2(:,k) = A * xupd_2(:,k-1);
%     Ppre_2(:,:,k) = A * Pupd_2(:,:,k-1) *A' + Bd*Q_cons(:,k-1)*Bd';
% 
%     if k<m
%         S_2(:,:,k)  = C*Ppre_2(:,:,k)*C' + diag(R_cons(:,k));
%         K_2(:,:,k) = Ppre_2(:,:,k)*C' / S_2(:,:,k);
%         R_hat2(:,k) = R_cons(:,k);
%     else % m+1
%         S_2(:,:,k)  = C*Ppre_2(:,:,k)*C' + R_hat_2(:,:,k);
%         K_2(:,:,k) = Ppre_2(:,:,k)*C' / S_2(:,:,k);
%     end
% 
%     xupd_2(:,k) = xpre_2(:,k) + K_2(:,:,k)*(Ytrue(:, k) - C*xpre_2(:,k) );
%     Pupd_2(:,:,k) = Ppre_2(:,:,k) - K_2(:,:,k)*C*Ppre_2(:,:,k); %
% 
%     e_y(:, k) =  Ytrue(:, k) -  C*xupd_2(:,k); %k=2才开始有值
% 
%     if k>=m
%         C_v = 0;
%         for j=1:m
%             C_v = C_v + e_y(:,k-j+1)*e_y(:,k-j+1)';
%         end
%         R_hat_2(:,:,k+1) = C_v/m + C*Pupd_2(:,:,k)*C'; % k时刻的估计给k+1时刻使用
%         R_hat2(:,k+1) = diag(R_hat_2(:,:,k+1),0);
%     end
%     RMS_2(1, k, mc) = sum((xupd_2(1:2, k) -Xtrue(1:2, k)).^2, 1);
%     RMS_2_1(1, k, mc) = sum((xupd_2(3:4, k) -Xtrue(3:4, k)).^2, 1);
% 
%     SRNFN_P2(1, k, mc) = 1/4^2*trace(( Pupd_2(:,:,k)  - Pupd(:,:,k) )*( Pupd_2(:,:,k)  - Pupd(:,:,k) )');
%     SRNFN_R2(1, k, mc) = 1/3^2*trace(( diag(R_hat2(:,k)) - diag(R_true(:, k)))*( diag(R_hat2(:,k)) - diag(R_true(:, k)))');
% 
% end
% 
% %% Stretch 3, E
% for k = 1002:1501 % A
%     xpre_2(:,k) = A * xupd_2(:,k-1);
%     Ppre_2(:,:,k) = A * Pupd_2(:,:,k-1) *A' + Bd*Q_cons(:,k-1)*Bd';
% 
%     if k<m
%         S_2(:,:,k)  = C*Ppre_2(:,:,k)*C' + diag(R_cons(:,k));
%         K_2(:,:,k) = Ppre_2(:,:,k)*C' / S_2(:,:,k);
%         R_hat2(:,k) = R_cons(:,k);
%     else % m+1
%         S_2(:,:,k)  = C*Ppre_2(:,:,k)*C' + R_hat_2(:,:,k);
%         K_2(:,:,k) = Ppre_2(:,:,k)*C' / S_2(:,:,k);
%     end
% 
%     xupd_2(:,k) = xpre_2(:,k) + K_2(:,:,k)*(Ytrue(:, k) - C*xpre_2(:,k) );
%     Pupd_2(:,:,k) = Ppre_2(:,:,k) - K_2(:,:,k)*C*Ppre_2(:,:,k); %
% 
%     e_y(:, k) =  Ytrue(:, k) -  C*xupd_2(:,k); %k=2才开始有值
% 
%     if k>=m
%         C_v = 0;
%         for j=1:m
%             C_v = C_v + e_y(:,k-j+1)*e_y(:,k-j+1)';
%         end
%         R_hat_2(:,:,k+1) = C_v/m + C*Pupd_2(:,:,k)*C'; % k时刻的估计给k+1时刻使用
%         R_hat2(:,k+1) = diag(R_hat_2(:,:,k+1),0);
%     end
%     RMS_2(1, k, mc) = sum((xupd_2(1:2, k) -Xtrue(1:2, k)).^2, 1);
%     RMS_2_1(1, k, mc) = sum((xupd_2(3:4, k) -Xtrue(3:4, k)).^2, 1);
% 
%     SRNFN_P2(1, k, mc) = 1/4^2*trace(( Pupd_2(:,:,k)  - Pupd(:,:,k) )*( Pupd_2(:,:,k)  - Pupd(:,:,k) )');
%     SRNFN_R2(1, k, mc) = 1/3^2*trace(( diag(R_hat2(:,k)) - diag(R_true(:, k)))*( diag(R_hat2(:,k)) - diag(R_true(:, k)))');
% 
% end