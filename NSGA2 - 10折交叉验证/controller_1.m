clc;clear;
rng default
warning('off')
%% parameter
ms=250; mu=30; ks=10000; kt=100000; cs=1000;

ts=0.01;%
T=10;   %
t = 0:ts:T;
%% plant
A = [0     0      1    -1;
    0     0      0     1;
    -ks/ms 0 -cs/ms cs/ms;
    ks/mu -kt/mu cs/mu -cs/mu];
Bu = [0 ;0 ;1/ms ;-1/mu];
Bd = [0 0;0 -1;0 0;0 0];

C = [-ks/ms 0 -cs/ms cs/ms;
    1 0 0 0;
    0 1 0 0];
Du = [1/ms;0;0];
Dd = zeros(3,2);
%离散化
sys = ss(A,[Bu Bd],C,[Du Dd],'inputname',{'u' 'w1' 'w2'},'outputname','y');
sys = c2d(sys,ts);

A = sys.A;B = sys.B;
Bu = B(:,1) ;Bd = sys.B(:,2:3);

Nx = size(A,1);Ny = size(C,1);Nu = size(Du,2);
%m、p
p = 20; %15
m = 3; %3
%% road
%road-1（single bump）
t0 = 2; %激励起点
width = 3; %激励周期/宽度
t1 = t0+width; %激励终点/因为只画一个波，所以只有一个周期内的数据
h = 5; %2h是高度/cm,计算时要/100，用m去计算
for i = 1:length(t)
    if t(i)>=t0 && t(i)<=t1
        w(i) = h*(1-cos((2*pi/width)*(t(i)-t0)))/100;
    else
        w(i) = 0;
    end  
end
% % Q = cov(w);
% 
% type = [16; 64; 256; 1024;4096;16384;65536;262144]*10^-6;
% speed=[10;20;30;40;50;60;70;80;90];
% % （A）
% fe=28.3;
% fs=0.33;
% N=500;
% 
% delta_f=(fe-fs)/N;
% sita = rand(N,1) * 2 * pi;
% %rand(m,1) * 2 * pi,rand是平均分布，不是标准正态分布
% for i=1:length(t)
%     w(i) = 0;
%     for j=1:N
%         f(j)=fs+delta_f*(j-0.5); %中心频率 
%         w(i)=w(i)+sqrt((2*type(1)*(0.1^2)*speed(3)*delta_f)/(f(j)^2))*sin(2*pi*f(j)*t(i)+sita(j,1));
%     end
% end
% w(1)=0;

%建立W、delta_W
for k = 1:length(t)-p
    %delta_w
    if k==1
%     if k==t0/ts+2 ||k==1
        delta_w(1)=0;
        for i=2:p
            delta_w(i) = w(k+i-1)-w(k+i-2); %delta_d(i,1)
        end
    else
        for i=1:p
            delta_w(i) = w(k+i-1)-w(k+i-2); %delta_d(i,1)
        end
    end
    %delta_w_dot
    if k==1 
%     if k==t0/ts+2 ||k==1
        delta_w_dot(1)=0;
        delta_w_dot(2)=0;
        for i=3:p
            delta_w_dot(i) = ( w(k+i-1)-2*w(k+i-2)+w(k+i-3) )/0.01; %delta_d(i,1)
        end
    elseif k==2 
%     elseif k==t0/ts+3 ||k==2
        delta_w_dot(1)=0;
        for i=2:p
            delta_w_dot(i) = ( w(k+i-1)-2*w(k+i-2)+w(k+i-3) )/0.01; %delta_d(i,1)
        end
    else
        for i=1:p
            delta_w_dot(i) = ( w(k+i-1)-2*w(k+i-2)+w(k+i-3) )/0.01; %delta_d(i,1)
        end
    end
    
    delta_W(:,k) = [delta_w';delta_w_dot'];
end
for k = 1:length(t)
   %w(k)
    W(1,k) =  w(k);
    %w_dot(k)
     if k==1
%       if k==1 %随机路面从k=1就有值
        W(2,k) = 0;
    else
        W(2,k) = ( w(k)-w(k-1) )/0.01;
    end  
end
%% mpc
%计算Sx,I,Sd,Su,Kmpc
A_cell = cell(p,1); %Sx
B1_cell = cell(p,m); %Su1
B2_cell = cell(p,m); %Su2
C1_cell = cell(p,p); %Sd1
C2_cell = cell(p,p); %Sd2
D_cell = cell(p,1);%I
for j = 1:1:p
    %Sx
    a1 = 0;
    for i = 1:1:j
        a1 = a1+ C * A^i;
    end %只有循环完才会跳出
    A_cell{j,1}=a1;
    
    %I
    D_cell{j,1}=eye(Ny);
    
    %Sd1
    for k=1:p %给每一列赋值
        c1 = 0;
        if k<=j
            for i=1:j-k+1 %叠加的符号
                c1 =c1+ C * (A^(i-1)) * Bd;
            end
        else
            c1 = zeros(size(C,1),size(Bd,2));%和Dd同阶
        end
        C1_cell{j,k}=c1;
    end
    %Sd2(p×p),应为书上是y=cx，这里是y=cx+du+dw
    for k=1:p %给每一列赋值
        c2 = 0;
        if k~=1&&k<=j+1
            c2 = c2 + Dd;
        else
            c2 = zeros(size(Dd,1),size(Dd,2));%和Dd同阶
        end
        C2_cell{j,k}=c2;
    end
    
    %Su1(p×m)
    for k=1:1:m %给每一列赋值
        b1 = 0;
        if k<=j
            for i=1:j-k+1
                b1 = b1+C * (A^(i-1)) * Bu;
            end
        else
            b1 = zeros(size(C,1),size(Bu,2));%和Cc*Bu同阶
        end
        B1_cell{j,k} = b1;
    end
    %Su2(p×m),应为书上是y=cx，这里是y=cx+du
    for k=1:m %给每一列赋值
        b2 = 0;
        if k~=1&&k<=j+1
            b2 = b2 + Du;
        else
            b2 = zeros(size(Du,1),size(Du,2));%和Dd同阶
        end
        B2_cell{j,k}=b2;
    end
end

Sx=cell2mat(A_cell);
Sd1=cell2mat(C1_cell);Sd2=cell2mat(C2_cell);Sd = Sd1+Sd2;
Su1=cell2mat(B1_cell);Su2=cell2mat(B2_cell);Su = Su1+Su2;
I =cell2mat(D_cell);

%中间量D1、D2
for i=1:m
    if i==1
        D1(1,i)=1;
    else
        D1(1,i)=0;
    end
end
for i=1:p
    if i==1
        D2(1,i)=1;
    else
        D2(1,i)=0;
    end
end
D2=blkdiag(D2,D2);

%  L_w=[1 50 500 5e-2];
L_w=[0.007120561	619.2157902	174.9943667	4005.397014];

Ly_i = [L_w(1) 0 0 ;0 L_w(2) 0;0 0 L_w(3)]; %Nc个输出变量
Lu_i = L_w(4); %Nu个控制变量

Ly=cell(p,p);
for i=1:1:p
    for j=1:1:p
        if i==j
            Ly{i,j} = Ly_i;
        else
            Ly{i,j} = zeros(Ny,Ny);
        end
    end
end
Ly=cell2mat(Ly);%权重矩阵

Lu=cell(m,m);
for i=1:1:m
    for j=1:1:m
        if i==j
            Lu{i,j} = Lu_i;
        else
            Lu{i,j} = zeros(Nu,Nu);
        end
    end
end
Lu=cell2mat(Lu);%权重矩阵

%参考
reff = cell(p,1);
for k=1:1:p
    reff{k,1}= zeros(3,1);
end
reff = cell2mat(reff);

%% 约束生成区域
%△U
A_T=cell(m,m);
for i=1:1:m
    for j=1:1:m
        if i==j
            A_T{i,j}=eye(Nu);
        else
            A_T{i,j}=zeros(Nu);
        end
    end
end
A_T = cell2mat(A_T);
delta_umin=-1000;
delta_umax=1000; 
delta_Umin=kron(ones(m,1),delta_umin); %给最大电流，阻尼器能阶跃的最大值
delta_Umax=kron(ones(m,1),delta_umax); %

%U
A_L = cell(m,m);% Nc=m
for i=1:1:m
    for j=1:1:m
        if j<=i
            A_L{i,j} = eye(Nu);
        else
            A_L{i,j} = 0;
        end
    end
end
A_L = cell2mat(A_L);
umin = -4000;%维数与控制变量的个数相同
umax = 4000;

%Y
y_max=[inf;0.1;0.15]; %
y_min=[inf;-0.1;-0.15]; %
y_max=kron(ones(p,1),y_max);
y_min=kron(ones(p,1),y_min);

%初始化
x(:,1) = [0;0;0;0];
xhat(:,1) = [0;0;0;0];

x2(:,1) = [0;0;0;0];
u2(:,1) = 0; 
y2(:,1) = C*x(:,1) + Dd*W(:,1) + Du*u2(:,1);

P(:,:) = zeros(4,4);
Q =  [ 5.062781109445248/10000, 0;
    0, 8.224369693589264/10000 ]; %w和w_dot两个参数的过程噪声
R = diag([0.01 , 0.01 , 0.01]).^2;

for k = 2:length(t)-p %新的时间0:Ts:T-p*Ts
    % passive（离散,并不是最优）
    x(:,k) = A*x(:,k-1) + Bd*W(:,k-1);
    y(:,k) = C*x(:,k) + Dd*W(:,k);
    %mpc
    x2(:,k) = A*x2(:,k-1) + Bd*W(:,k-1) + Bu*u2(k-1) ;
  
    %% mpc
    delta_x = x2(:,k)- x2(:,k-1);

    %U
    Umin = kron(ones(m,1), -umin + u2(k-1)); 
    Umax = kron(ones(m,1), -u2(k-1) + umax); 
    %Y
    U_Ymax = -( (Sx+I*C)*delta_x + I*y2(:,k-1)  + (Sd+I*Dd*D2)*delta_W(:,k) ) + y_max;
    U_Ymin = ( (Sx+I*C)*delta_x + I*y2(:,k-1)  + (Sd+I*Dd*D2)*delta_W(:,k) ) - y_min;
  
    Su_max = Su + I*Du*D1; %1×m
    for l=1:size(Su,1)
        if mod(l,3)==1
            Su_max(l,:) = 0;
            U_Ymax(l,:) = 0;
        end
    end
    Su_min = Su + I*Du*D1; %1×m
    for h=1:size(Su,1)
        if mod(h,3)==1
            Su_min(h,:) = 0;
            U_Ymin(h,:) = 0;
%         elseif mod(h,3)==0
%             Su_min(h,:) = 0;
%             U_Ymin(h,:) = 0;
        end
    end
    %% 结合控制量约束和输出量约束
%     a = [A_T; -A_T]; %1
%     b = [delta_Umax; delta_Umin];
%     
    a = [A_L; -A_L]; %2
    b = [ Umax; Umin];
%     
%     a = [ Su_max; -Su_min]; %3
%     b = [ U_Ymax; U_Ymin];
%     
%     a = [A_T; -A_T; A_L; -A_L]; %1,2
%     b = [delta_Umax; delta_Umin; Umax; Umin];
%     
%     a = [A_T; -A_T; Su_max; -Su_min]; %1,3
%     b = [delta_Umax; delta_Umin; U_Ymax; U_Ymin];
%     
%      a = [A_L; -A_L; Su_max; -Su_min]; %2,3
%      b = [Umax; Umin; U_Ymax; U_Ymin];
%     
%     a = [A_T; -A_T; A_L; -A_L; Su_max; -Su_min]; %1,2,3
%     b = [delta_Umax; delta_Umin; Umax; Umin; U_Ymax; U_Ymin];

    %%  求解
    H = 2*( (Su + I*Du*D1)' *(Ly'*Ly) * (Su + I*Du*D1) + (Lu' * Lu));
    G = -2*( (Su + I*Du*D1)' *(Ly'*Ly) * (reff - (Sx+I*C)*delta_x - I*y2(:,k-1) - (Sd+I*Dd*D2)*delta_W(:,k)));
    %二次规划中△u的初始化
    for i=1:m
        x0(i,1) = rand();
    end
    
    options = optimoptions('quadprog','Algorithm','interior-point-convex','Display','off');
%      options = optimoptions('quadprog','Algorithm','interior-point-convex');
%       U = quadprog(H, G, a, b,[],[],[],[],x0,options);
     U = quadprog(H, G, a, b,[],[],[],[],[],options);
%     U = quadprog(H, G, [], [],[],[],[],[],[],options);
    delta_u2 = U(1);
    sanjiao_u2(1,k) = delta_u2;
    u2(:,k) = u2(:,k-1) + delta_u2;
    y2(:,k) = C*x2(:,k) + Du*u2(:,k) + Dd*W(:,k);

end
t = 0:ts:T-p*ts;
figure
plot(t,u2)
title('L_w=[0.1 200 1 1e-2]')
figure
plot(t,y2(1,:),t,y(1,:));legend('mpc','passive')
figure
plot(t,y2(2,:),t,y(2,:));legend('mpc','passive')
figure
plot(t,y2(3,:),t,y(3,:));legend('mpc','passive')
