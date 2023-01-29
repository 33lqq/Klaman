clc;clear
% function nsga_2_optimization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  p=parpool(20);
tic
pop = 150; %种群数量
gen = 200; %迭代次数
M = 4; %目标函数数量，f1、f2
V = 226; %维度,Xi的变量个数，f(x) = 
% x1 + x2 + ... + xV
%268,52
% range = 2e3*ones(1,V);

min_range = zeros(1,V); %下界,生成1*100的个体向量，全为0
max_range = ones(1,V); %上界,生成1*100的个体向量，全为1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%所有样本特征参数矩阵合成一个矩阵
load('A.mat');load('G.mat');load('H.mat');load('F.mat');load('E.mat');load('C.mat');load('B.mat');load('D.mat')
% load('A_1.mat');load('G_1.mat');load('H_1.mat');load('F_1.mat');load('E_1.mat');load('C_1.mat');load('B_1.mat');load('D_1.mat')
% load('A_2.mat');load('G_2.mat');load('H_2.mat');load('F_2.mat');load('E_2.mat');load('C_2.mat');load('B_2.mat');load('D_2.mat')

data_num = 200;
level = ['A','B','C','D','E','F','G','H'];
for i=1:8
%     data_1 = eval(level(i));
    data(data_num*(i-1)+1:data_num*i,:) = eval(level(i));
    label(data_num*(i-1)+1:data_num*i,:) = ones(data_num,1)*i;
end

%归一化
[data,inputps]=mapminmax(data);
randIndex = randperm(size(data,1));
data = data(randIndex,:);
label = label(randIndex,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
chromosome = initialize_variables(pop, M, V, min_range, max_range, data, label); %初始化种群
% chromosome是一个矩阵，行数是pop的个数，列数是(V+M+2=决策变量数+目标函数维数+等级+拥挤度)
chromosome = non_domination_sort_mod(chromosome, M, V); %对初始化种群进行非支配快速排序和拥挤度计算

f1_opt = zeros(1,gen);
f2_opt = zeros(1,gen);
f3_opt = zeros(1,gen);
f4_opt = zeros(1,gen);
for i = 1 : gen
    pool = round(pop/2); %round() 四舍五入取整 交配池大小 ,100
    tour = 2; %竞标赛  参赛选手个数
    parent_chromosome = tournament_selection(chromosome, pool, tour); %竞标赛选择适合繁殖的父代
    mu = 20; %交叉和变异算法的分布指数
    mum = 20;
    %进行交叉变异产生子代
    offspring_chromosome = genetic_operator(parent_chromosome,M, V, mu, mum, min_range, max_range, data, label);
    [main_pop,~] = size(chromosome); %父代种群的大小
    [offspring_pop,~] = size(offspring_chromosome); %子代种群的大小
    
    clear temp
    intermediate_chromosome(1:main_pop,:) = chromosome;
    %合并父代和子代种群
    intermediate_chromosome(main_pop + 1 : main_pop + offspring_pop,1 : M+V) = offspring_chromosome;
    %对新的种群进行非支配排序
    intermediate_chromosome = non_domination_sort_mod(intermediate_chromosome, M, V);
   %选择排序前N个优秀的个体构成新种群
    chromosome = replace_chromosome(intermediate_chromosome, M, V, pop);
    
    %% 统计每轮中最优的适应度值
    index = find(chromosome(:,V+6)==inf); %是否在等级1中统计？,inf是每个变量的最小值
    chromo =  chromosome(intersect(1:pop, index),:); %union并集，intersect交集
    
    f1_opt(:,i) = min(chromo(:,V+1));
    f2_opt(:,i) = min(chromo(:,V+2));
    f3_opt(:,i) = min(chromo(:,V+3));
    f4_opt(:,i) = min(chromo(:,V+4));
      
    if ~mod(i,100)
        clc;
        fprintf('%d generations completed\n',i);
    end
    i
end

%% plot
figure
subplot(2,2,1)
plot(1:gen,f1_opt);
subplot(2,2,2)
plot(1:gen,f2_opt);
subplot(2,2,3)
plot(1:gen,f3_opt);
subplot(2,2,4)
plot(1:gen,f4_opt);

figure
if M == 2
    plot(chromosome(:,V + 1),chromosome(:,V + 2),'*');
    xlabel('f_1'); ylabel('f_2');
    title('Pareto Optimal Front');
elseif M == 3
    plot3(chromosome(:,V + 1),chromosome(:,V + 2),chromosome(:,V + 3),'*');
%     ax = gca;
%     ax.XDir = 'reverse'; 
%     view(-31,14) % 相机视角
    xlabel('f_1'); ylabel('f_2'); zlabel('f_3');
    title('Pareto Optimal Surface');
elseif M == 4
    figure(1)
    set(gcf,'unit','centimeters','position',[10,10,12,7.5])    % 图形窗口在电脑屏幕上的位置和尺寸[左 下 宽 高]
    linewidth_line = 1.5;      % 图形线条宽度
    markersize = 4;            % 图形标记点大小
    linewidth_gca = 0.7;       % 横纵坐标轴宽度
    fontsize_gca = 10;         % 横纵坐标轴刻度字体大小
    fontsize_label = 12;       % 横纵坐标轴字体大小
   
    f1 = chromosome(:,V+1);
    f3 = chromosome(:,V+3);
    f4 = chromosome(:,V+4);
    f2 = chromosome(:,V+2);

%     f1 = chromosome(:,V+1);
%     f3 = chromosome(:,V+3);
%     f4 = chromosome(:,V+4);
%     f2 = chromosome(:,V+2);
    scatter3(0.4+f1,0.005*f2,f3,40,f2,"filled") %40 圆点大小
%     ax = gca;
%     ax.XDir = 'reverse'; 
    view(-31,14) % 相机视角
    
    set(gca,'linewidth',linewidth_gca,'fontsize',fontsize_gca)
    set(gca,'GridLineStyle','--');
    xlabel('Classification model accuracy','fontsize',fontsize_label);
    ylabel('Feature inter-cluster distance','fontsize',fontsize_label);
    zlabel('Feature intra-cluster distance','fontsize',fontsize_label);
    
    h = rotate3d;
    set(h, 'ActionPreCallback', 'set(gcf,''windowbuttonmotionfcn'',@align_axislabel)')
    set(h, 'ActionPostCallback', 'set(gcf,''windowbuttonmotionfcn'','''')')
    set(gcf, 'ResizeFcn', @align_axislabel)
    align_axislabel([], gca)
    axislabel_translation_slider;
    
    cb = colorbar;
    cb.Label.String = 'Number of feature subsets';
    
    hfig = figure(1);
    figWidth = 12;  % 设置图片宽度
    figHeight = 7.5;  % 设置图片高度
    set(hfig,'PaperUnits','centimeters'); % 图片尺寸所用单位
    set(hfig,'PaperPosition',[0 0 figWidth figHeight]);
    fileout = 'naga2.'; % 输出图片的文件名
    % print(hfig,[fileout,'tif'],'-r300','-dtiff'); % 设置图片格式、分辨率
    print(hfig,[fileout,'pdf'],'-r600','-dpdf'); % 设置图片格式、分辨率
end
toc
% delete(p);