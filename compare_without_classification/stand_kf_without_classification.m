% coede for:
% 标准卡尔曼
clc;clear;close all;
dbstop if error;
% seed = 1;
% rng(seed)

simTime = 15;                       % [s]
Ts = 0.01; %采样时间，这点和车速挂钩，
time=0:Ts:simTime;
N=length(time);% 25s+1

%% plant
ms=250; mu=30; ks=10000; kt=100000; cs=1000;

A = [0 0 1 -1;
    0 0 0 1;
    -ks/ms 0 -cs/ms cs/ms;
    ks/mu -kt/mu cs/mu -cs/mu];
Bd = [0 -1 0 0]';
C = [-ks/ms 0 -cs/ms cs/ms;
    1 0 0 0;
    0 1 0 0];

sys = ss(A,Bd,C,[]);sys = c2d(sys,Ts);

A = sys.A;
Bd = sys.B; %两种方法A基本不变，Bd会改变 %B不是输入u，是不是不需要离散化
C = sys.C;

%% Time_vary
[Q,R,Q_cons,Q_true,R_cons,R_true] = PNCM_RNCM(simTime, Ts);

%% Simulation
mcTime = 2000;
RMS = zeros(1, N, mcTime);
RMS1 = zeros(1, N, mcTime);
RMS_1 = zeros(1, N, mcTime); RMS_1_1 = zeros(1, N, mcTime); SRNFN_P1 = zeros(1, N, mcTime); SRNFN_R1 = zeros(1, N, mcTime);
RMS_2 = zeros(1, N, mcTime); RMS_2_1 = zeros(1, N, mcTime); SRNFN_P2 = zeros(1, N, mcTime); SRNFN_R2 = zeros(1, N, mcTime);
RMS_3 = zeros(1, N, mcTime); RMS_3_1 = zeros(1, N, mcTime); SRNFN_P3 = zeros(1, N, mcTime); SRNFN_R3 = zeros(1, N, mcTime);
RMS_4 = zeros(1, N, mcTime); RMS_4_1 = zeros(1, N, mcTime); SRNFN_P4 = zeros(1, N, mcTime); SRNFN_R4 = zeros(1, N, mcTime); 
RMS_41 = zeros(1, N, mcTime); RMS_41_1 = zeros(1, N, mcTime); SRNFN_P41 = zeros(1, N, mcTime); SRNFN_R41 = zeros(1, N, mcTime);

for mc= 1:mcTime
    mc
    %% road && measurement noise
    for i=1:N
        %         w(i,3) = mvnrnd(zeros(1, 1), Q_true(:,i))';
        w(i,1) = i*0.01;
        w(i,2) = mvnrnd(zeros(1, 1), Q_true(:,i))';
    end
    for i=1:N
        v(:,i) = mvnrnd(zeros(3, 1), diag(R_true(:,i)))';
    end
    %% 真实值
    Xtrue = zeros(4, N); %用来对比观测的状态变量
    Xtrue(:, 1) = zeros(4,1); %Xtrue（0）列向量
    Ytrue = zeros(3, N);
    for k = 2:N %A路面
        Xtrue(:, k) = A * Xtrue(:, k-1) + Bd * w(k-1,2);
    end

    for k = 1:N
        Ytrue(:, k) = C*Xtrue(:, k) + v(:,k);
    end
    %% 真实Q，R
    [xupd, Pupd, RMS, RMS1] = KFTCM_without_classification(A,Bd,C,Q_true,R_true,Xtrue,Ytrue,mc, RMS, RMS1);
    %% 固定Q，R
    [xupd_1, Pupd_1, RMS_1, RMS_1_1, SRNFN_P1, SRNFN_R1] = KFNCM_without_classification(A,Bd,C,Q_cons,R_cons,R_true,Xtrue,Ytrue,Pupd,mc, RMS_1, RMS_1_1, SRNFN_P1, SRNFN_R1);
    %% 方差匹配
    [xupd_2, Pupd_2, R_hat2, RMS_2, RMS_2_1, SRNFN_P2, SRNFN_R2] = RSAKF_without_classification(A,Bd,C,Q_cons,R_cons,R_true,Xtrue,Ytrue,Pupd,mc, RMS_2, RMS_2_1, SRNFN_P2, SRNFN_R2);
    %% VB-AKF-R
    [xupd_3, Pupd_3, R_hat3, RMS_3, RMS_3_1, SRNFN_P3, SRNFN_R3] = VBAKF_R_without_classification(A,Bd,C,Q,R,Q_cons,R_cons,R_true,Xtrue,Ytrue,Pupd,mc,  RMS_3, RMS_3_1, SRNFN_P3, SRNFN_R3);
    %% 变分贝叶斯
    [xupd_4, Pupd_4, R_hat4, R_hat41, RMS_4, RMS_4_1, RMS_41, RMS_41_1, SRNFN_P4, SRNFN_R4, SRNFN_P41, SRNFN_R41] = My_filter_without_classification(A,Bd,C,Q,R,Q_cons,R_cons,R_true,Xtrue,Ytrue,Pupd,mc, RMS_4, RMS_4_1, RMS_41, RMS_41_1, SRNFN_P4, SRNFN_P41, SRNFN_R4,SRNFN_R41);
end



