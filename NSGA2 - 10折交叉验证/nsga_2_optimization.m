clc;clear
% function nsga_2_optimization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  p=parpool(20);
tic
pop = 150; %��Ⱥ����
gen = 200; %��������
M = 4; %Ŀ�꺯��������f1��f2
V = 226; %ά��,Xi�ı���������f(x) = 
% x1 + x2 + ... + xV
%268,52
% range = 2e3*ones(1,V);

min_range = zeros(1,V); %�½�,����1*100�ĸ���������ȫΪ0
max_range = ones(1,V); %�Ͻ�,����1*100�ĸ���������ȫΪ1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%��������������������ϳ�һ������
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

%��һ��
[data,inputps]=mapminmax(data);
randIndex = randperm(size(data,1));
data = data(randIndex,:);
label = label(randIndex,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
chromosome = initialize_variables(pop, M, V, min_range, max_range, data, label); %��ʼ����Ⱥ
% chromosome��һ������������pop�ĸ�����������(V+M+2=���߱�����+Ŀ�꺯��ά��+�ȼ�+ӵ����)
chromosome = non_domination_sort_mod(chromosome, M, V); %�Գ�ʼ����Ⱥ���з�֧����������ӵ���ȼ���

f1_opt = zeros(1,gen);
f2_opt = zeros(1,gen);
f3_opt = zeros(1,gen);
f4_opt = zeros(1,gen);
for i = 1 : gen
    pool = round(pop/2); %round() ��������ȡ�� ����ش�С ,100
    tour = 2; %������  ����ѡ�ָ���
    parent_chromosome = tournament_selection(chromosome, pool, tour); %������ѡ���ʺϷ�ֳ�ĸ���
    mu = 20; %����ͱ����㷨�ķֲ�ָ��
    mum = 20;
    %���н����������Ӵ�
    offspring_chromosome = genetic_operator(parent_chromosome,M, V, mu, mum, min_range, max_range, data, label);
    [main_pop,~] = size(chromosome); %������Ⱥ�Ĵ�С
    [offspring_pop,~] = size(offspring_chromosome); %�Ӵ���Ⱥ�Ĵ�С
    
    clear temp
    intermediate_chromosome(1:main_pop,:) = chromosome;
    %�ϲ��������Ӵ���Ⱥ
    intermediate_chromosome(main_pop + 1 : main_pop + offspring_pop,1 : M+V) = offspring_chromosome;
    %���µ���Ⱥ���з�֧������
    intermediate_chromosome = non_domination_sort_mod(intermediate_chromosome, M, V);
   %ѡ������ǰN������ĸ��幹������Ⱥ
    chromosome = replace_chromosome(intermediate_chromosome, M, V, pop);
    
    %% ͳ��ÿ�������ŵ���Ӧ��ֵ
    index = find(chromosome(:,V+6)==inf); %�Ƿ��ڵȼ�1��ͳ�ƣ�,inf��ÿ����������Сֵ
    chromo =  chromosome(intersect(1:pop, index),:); %union������intersect����
    
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
%     view(-31,14) % ����ӽ�
    xlabel('f_1'); ylabel('f_2'); zlabel('f_3');
    title('Pareto Optimal Surface');
elseif M == 4
    figure(1)
    set(gcf,'unit','centimeters','position',[10,10,12,7.5])    % ͼ�δ����ڵ�����Ļ�ϵ�λ�úͳߴ�[�� �� �� ��]
    linewidth_line = 1.5;      % ͼ���������
    markersize = 4;            % ͼ�α�ǵ��С
    linewidth_gca = 0.7;       % ������������
    fontsize_gca = 10;         % ����������̶������С
    fontsize_label = 12;       % ���������������С
   
    f1 = chromosome(:,V+1);
    f3 = chromosome(:,V+3);
    f4 = chromosome(:,V+4);
    f2 = chromosome(:,V+2);

%     f1 = chromosome(:,V+1);
%     f3 = chromosome(:,V+3);
%     f4 = chromosome(:,V+4);
%     f2 = chromosome(:,V+2);
    scatter3(0.4+f1,0.005*f2,f3,40,f2,"filled") %40 Բ���С
%     ax = gca;
%     ax.XDir = 'reverse'; 
    view(-31,14) % ����ӽ�
    
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
    figWidth = 12;  % ����ͼƬ���
    figHeight = 7.5;  % ����ͼƬ�߶�
    set(hfig,'PaperUnits','centimeters'); % ͼƬ�ߴ����õ�λ
    set(hfig,'PaperPosition',[0 0 figWidth figHeight]);
    fileout = 'naga2.'; % ���ͼƬ���ļ���
    % print(hfig,[fileout,'tif'],'-r300','-dtiff'); % ����ͼƬ��ʽ���ֱ���
    print(hfig,[fileout,'pdf'],'-r600','-dpdf'); % ����ͼƬ��ʽ���ֱ���
end
toc
% delete(p);