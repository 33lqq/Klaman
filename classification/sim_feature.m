%% 生成特征数据
% coede for:
% 标准卡尔曼
clc;clear;close all;
dbstop if error;
% seed = 1;
% rng(seed)

simTime = 15;
n = 4;
Ts = 0.01; %采样时间，这点和车速挂钩，
t=0:Ts:simTime;
N=length(t);% 25s+1

%% plant
ms=250; mu=30; ks=10000; kt=100000; cs=1000;

A = [0 0 1 -1;
    0 0 0 1;
    -ks/ms 0 -cs/ms cs/ms;
    ks/mu -kt/mu cs/mu -cs/mu];
% B = [0 0 0 0]';
Bd = [0 -1 0 0]';
C = [-ks/ms 0 -cs/ms cs/ms;
    1 0 0 0;
    0 1 0 0];


sys = ss(A,Bd,C,[]);sys = c2d(sys,Ts);

A = sys.A;
Bd = sys.B; %两种方法A基本不变，Bd会改变 %B不是输入u，是不是不需要离散化
C = sys.C;

number = 8;
level = ['A','B','C','D','E','F','G','H'];
%% Time_vary
a = 1; % 超调
b = 4; %周期

Q = [0.0022,0.0088,0.0351,0.1406,0.5624,2.2495,8.9979,35.992]; %过程噪声协方差矩阵,6×1
% Q = ones(4,1)*Q;
Q1 = Q(:,number)*ones(1,N);%A
for k=1:N
    Q_1(:,k) = (1 - a*sin(b*pi*k/N))*Q1(:,k);
end
Q_true = Q_1;

R =[0.0001,0.0001,0.0015;0.0004,0.0004,0.006;0.0016,0.0016,0.024;...
    0.0064,0.0064,0.096;0.0256,0.0256,0.384;0.1024,0.1024,1.536;...
    0.4096,0.4096,6.144;1.6384,1.6384,24.576]';
R1 = R(:,number)*ones(1,N);%A
for k=1:N
    for i=1:3
        if i==2
            R_1(i,k) = (1 - a*sin(b*pi*k/N))*R1(i,k);
        else
            R_1(i,k) = (1 + a*sin(b*pi*k/N))*R1(i,k);
        end
    end
end
R_true = [R_1];

mcTime = 10000;
RMS = zeros(1, N, mcTime);

speed = (40:40);

file_path = ['D:/MATLAB/kalman/feature/',level(number)];
if ~isfolder(path)
    mkdir(file_path); %文件夹是否存在，不存在则创建
end
s = dir(file_path);
if length(s) == 2 %文件为空
    disp("文件夹为空")
else
    delete(['D:/MATLAB/kalman/feature/',level(number),'/*.csv'])
    disp("文件夹已清空")
end
% for speed = (datasample(speed, 1,'Replace',false))

for mc= 1:mcTime
    mc
    %% 噪声信号
    for i=1:N
        w(:,i) = mvnrnd(zeros(1, 1), Q_true(:,i))';
        v(:,i) = mvnrnd(zeros(3, 1), diag(R_true(:,i)))';
    end
    %% 真实值
    Xtrue = zeros(4, N); %用来对比观测的状态变量
    Xtrue(:, 1) = zeros(4,1); %Xtrue（0）列向量
    Ytrue = zeros(3, N);
    for k = 2:N %A路面
        Xtrue(:, k) = A * Xtrue(:, k-1) + Bd * w(:,k-1);
        %         Xtrue(:, k) = A * Xtrue(:, k-1) + Bd * mvnrnd(zeros(1, 1), Q_cons(:,k-1))';
    end

    for k = 1:N
        Ytrue(:, k) = C*Xtrue(:, k) + v(:,k);
        %         Ytrue(:, k) = C*Xtrue(:, k) + mvnrnd(zeros(3, 1), diag(R_cons(:,k)))';
    end
    Ytrue(:,1)=[];
    Ytrue = Ytrue';

    writematrix(Ytrue, ['D:/MATLAB/kalman/feature/',level(number),'/',level(number),'_' ,num2str(speed),'_', num2str(mc,'%03d'),'.csv']);     % 写入csv

end

disp('数据生成完毕')