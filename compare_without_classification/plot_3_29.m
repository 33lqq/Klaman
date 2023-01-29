figure
hold on

% t0= 1;
% tf= 501;
% t0= 502;
% tf= 1001;
t0= 1002;
tf= 1501;

index = 3;
% plot((1:N),R_1(2, 1:N), '-r.')
plot((t0:tf),R_true(index, t0:tf))
plot((t0:tf),R_hat2(index, t0:tf))
plot((t0:tf),R_hat3(index, t0:tf))
plot((t0:tf),R_hat4(index, t0:tf))
% plot((t0:tf),R_hat41(2, t0:tf))
% plot((1:N),R_hat_s1(2, 1:N), '-g.')
legend('True', 'RSAKF','VBAKF_R', 'My_filter','Location', 'west');
% 
rmse1 = sqrt(mean(( R_hat1(2, 1:499)-R_1(2, 1:499) ).^2))
% rmse1_1 = sqrt(mean(( R_hat_s1(2, 1:499)-R_1(2, 1:499) ).^2))
% 
% SRNFN1 = ((1/3*mean(( R_hat1(2, 1:499)-R_1(2, 1:499) ).^2))).^(1/4)
% SRNFN1_1 = ((1/3*mean(( R_hat_s1(2, 1:499)-R_1(2, 1:499) )).^2)).^(1/4)
% 
% rmse1 = sqrt(mean(( xupd(2, 1:499)-Xtrue(2, 1:499) ).^2))
% rmse1_1 = sqrt(mean(( xupd_bwd(2, 1:499) - Xtrue(2, 1:499) ).^2))
% 
% 
% 
% for i=1:5
%     if i<3
%          i
%     else
%         break
%     end
% end

%% x对比


figure
hold on;
grid on;
plot(t(1:N), xupd(2, 1:N), '-r.')
plot(t(1:N), xupd_1(2, 1:N), '-B.')
plot(t(1:N), xupd_3(2, 1:N), '-Y.')
plot(t(1:N), xupd_4(2, 1:N))
plot(t(1:N), Xtrue(2, 1:N), '-k.')
legend('时变','固定','变分','实际')
xlabel('time-[s]')
ylabel('Tyre deflection-[m]')
title('轮胎行程')



figure
hold on;
grid on;
plot(t(1:N), xupd(3, 1:N), '-r.')
plot(t(1:N), xupd_1(3, 1:N), '-B.')
plot(t(1:N), xupd_3(3, 1:N), '-Y.')
plot(t(1:N), Xtrue(3, 1:N), '-k.')
legend('时变','固定','变分','实际')
xlabel('time-[s]')
ylabel('Sprung mass velocity-[m]')
title('簧载质量速度')

figure
hold on;
grid on;
plot(t(1:N), xupd(4, 1:N))
plot(t(1:N), xupd_1(4, 1:N))
plot(t(1:N), xupd_3(4, 1:N))
plot(t(1:N), xupd_4(4, 1:N))
plot(t(1:N), Xtrue(4, 1:N))
legend('时变','固定','变分','VB-R','实际')
xlabel('time-[s]')
ylabel('Unsprung mass velocity-[m]')
title('非簧载质量速度')

figure
hold on;
grid on;
plot(t(200:500), xupd(2, 200:500), '-r.')
% plot(t(200:500), xupd_bwd(2, 200:500), '-b.')
plot(t(200:500), Xtrue(2, 200:500), '-k.')
figure
hold on;
grid on;
plot(t(200:500), xupd(3, 200:500), '-r.')
% plot(t(200:500), xupd_bwd(3, 200:500), '-b.')
plot(t(200:500), Xtrue(3, 200:500), '-k.')
figure
hold on;
grid on;
plot(t(200:500), xupd(4, 200:500), '-r.')
% plot(t(200:500), xupd_bwd(4, 200:500), '-b.')
plot(t(200:500), Xtrue(4, 200:500), '-k.')
%% y对比
figure
hold on;
grid on;
plot(t(200:500), yupd(1, 200:500), '-r.')
% plot(t(200:500), xupd_bwd(1, 200:500), '-b.')
plot(t(200:500), Ytrue(1, 200:500), '-k.')
figure
hold on;
grid on;
plot(t(200:500), yupd(2, 200:500), '-r.')
% plot(t(200:500), xupd_bwd(1, 200:500), '-b.')
plot(t(200:500), Ytrue(2, 200:500), '-k.')
figure
hold on;
grid on;
plot(t(200:500), yupd(3, 200:500), '-r.')
% plot(t(200:500), xupd_bwd(1, 200:500), '-b.')
plot(t(200:500), Ytrue(3, 200:500), '-k.')

fprintf('滤波后误差均值x1为：%5.3f；均方差为：%5.3f\n', ...
    mean(abs(xupd(1,200:500)-Xtrue(1,200:500))), std(xupd(1,200:500)-Xtrue(1,200:500)));
fprintf('滤波后误差均值x2为：%5.3f；均方差为：%5.3f\n', ...
    mean(abs(xupd(2,200:500)-Xtrue(2,200:500))), std(xupd(2,200:500)-Xtrue(2,200:500)));
fprintf('滤波后误差均值x3为：%5.3f；均方差为：%5.3f\n', ...
    mean(abs(xupd(3,200:500)-Xtrue(3,200:500))), std(xupd(3,200:500)-Xtrue(3,200:500)));
fprintf('滤波后误差均值x4为：%5.3f；均方差为：%5.3f\n', ...
    mean(abs(xupd(4,200:500)-Xtrue(4,200:500))), std(xupd(4,200:500)-Xtrue(4,200:500)));

clc;clear
A = [0 1;-5 -2];
B = [0;3];
C = [0 1];
D = 0;

Ts = 0.25;
sys = ss(A,B,C,D);
sys1 = c2d(sys,Ts,'imp');  %A不变，B会改变，并且D也改变了
sys2 = c2d(sys,Ts); 
% sys = ss(A,B,C,D,Ts); %和连续系统输出一样
% A_1 = sys.A
% B_1 = sys.B %两种方法A基本不变，Bd会改变
% C_1 = sys.C
x(:,1) =[1 1]'; %0
% for i=2:10 %0.25
for i=2:10 %0.25
    x(:,i) = sys.A*x(:,i-1) ;
    y(:,i) = sys.C * x(:,i);
end

% Physical parameters
mb = 300;    % kg
mw = 60;     % kg
bs = 1000;   % N/m/s
ks = 16000 ; % N/m
kt = 190000; % N/m
% mb = 250;    % kg
% mw = 30;     % kg
% bs = 1000;   % N/m/s
% ks = 10000 ; % N/m
% kt = 100000; % N/m

% State matrices
A = [ 0 1 0 0; [-ks -bs ks bs]/mb ; ...
      0 0 0 1; [ks bs -ks-kt -bs]/mw];
B = [ 0 0; 0 1e3/mb ; 0 0 ; [kt -1e3]/mw];
C = [1 0 0 0; 1 0 -1 0; A(2,:)];
D = [0 0; 0 0; B(2,:)];

qcar = ss(A,B,C,D);
Ts = 0.25;
sys = c2d(qcar,Ts); 

figure
plot(Q_true,'*')
hold on
plot(Q_cons)
figure
plot(R_true(1,:),'*')
hold on
plot(R_cons(1,:))

figure
plot(w(1:499,2))
hold on
plot(w(1:499,3))
legend('谐波','正太')

%% 方差匹配
    m=20;
    xupd_2(:,1) = zeros(4,1);
    Pupd_2(:,:,1) = zeros(4);%eye(4);%zeros(4); %对称就行
    R_hat(:,:,m+1) = diag(R_cons(:,2));
    for k = 2:N %epoch k
        xpre_2(:,k) = A * xupd_2(:,k-1);
        Ppre_2(:,:,k) = A * Pupd_2(:,:,k-1) *A' + Bd*Q_cons(:,k-1)*Bd';

        if k<=m
            S_2(:,:,k)  = C*Ppre_2(:,:,k)*C' + diag(R_cons(:,k));
            K_2(:,:,k) = Ppre_2(:,:,k)*C' / S_2(:,:,k);
        else % m+1
            S_2(:,:,k)  = C*Ppre_2(:,:,k)*C' + R_hat(:,:,k);
            K_2(:,:,k) = Ppre_2(:,:,k)*C' / S_2(:,:,k);
        end

        xupd_2(:,k) = xpre_2(:,k) + K_2(:,:,k)*(Ytrue(:, k) - C*xpre_2(:,k) );
        Pupd_2(:,:,k) = Ppre_2(:,:,k) - K_2(:,:,k)*C*Ppre_2(:,:,k); %

        e_y(:, k) =  Ytrue(:, k) -  C*xupd_2(:,k); %k=2才开始有值

        if k<=m
            C_v = 0;
            for j=1:k
                C_v = C_v + e_y(:,k-j)*e_y(:,k-j)';
            end
            R_hat(:,:,k+1) = C_v/m + C*Pupd_2(:,:,k)*C'; % k时刻的估计给k+1时刻使用
        end
        if k>m
            C_v = 0;
            for j=1:m
                C_v = C_v + e_y(:,k-j)*e_y(:,k-j)';
            end
            R_hat(:,:,k+1) = C_v/m + C*Pupd_2(:,:,k)*C'; % k时刻的估计给k+1时刻使用
        end
        RMS_2(1, k, mc) = sum((xupd_2(1:2, k) -Xtrue(1:2, k)).^2, 1);
    end