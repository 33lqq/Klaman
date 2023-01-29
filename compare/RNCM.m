function [R_cons,R_true] = RNCM(simTime, Ts)
a = 1; % 超调
b = 4; %周期
n = 3;
L = fix(simTime/n);

R =[0.0001,0.0001,0.0015;0.0004,0.0004,0.006;0.0016,0.0016,0.024;...
    0.0064,0.0064,0.096;0.0256,0.0256,0.384;0.1024,0.1024,1.536;...
    0.4096,0.4096,6.144;1.6384,1.6384,24.576]';
% Stretch 1, A
x_1  = 0:Ts:L ;                 % 5[s]
R1 = R(:,1)*ones(1,length(x_1));%A
for k=1:length(x_1)
    R_1(:,k) = (1 + a*sin(b*pi*k/length(x_1)))*R1(:,k);
end
% Stretch 2, C
x_2  = L:Ts:2*L ;  
R2 = R(:,3)*ones(1,length(x_2));%C
for k=1:length(x_2)
    R_2(:,k) = (1 + a*sin(b*pi*k/length(x_2)))*R2(:,k);
end
% Stretch 3, E
x_3  = 2*L:Ts:3*L ; 
R3 = R(:,5)*ones(1,length(x_3));%E
for k=1:length(x_3)
    R_3(:,k) = (1 + a*sin(b*pi*k/length(x_3)))*R3(:,k);
end
R_cons = [R1,R2(:,2:end),R3(:,2:end)];
R_true = [R_1,R_2(:,2:end),R_3(:,2:end)];

%% plot(1)
figure(1)
hold on; box on; grid on;
set(gcf,'unit','centimeters','position',[10,10,8,5])    % 图形窗口在电脑屏幕上的位置和尺寸[左 下 宽 高]
% linewidth_line = 1.2;      % 图形线条宽度
linewidth_gca = 0.7;      % 横纵坐标轴宽度
fontsize_gca = 7;           % 横纵坐标轴刻度字体大小
fontsize_label = 9;         % 横纵坐标轴字体大小
fontsize_legend = 7;      % 图例字体大小

plot(R_cons(1,:),'LineWidth',2)
plot(R_true(1,:),'LineWidth',2)


line([501 501],[0 1.2],'linestyle','--','Color','k');
line([1001 1001],[0 1.2],'linestyle','--','Color','k');
text(200-10,0.9*1.2,'ISO-A','FontSize',fontsize_legend);
text(700-10,0.9*1.2,'ISO-C','FontSize',fontsize_legend);
text(1200-10,0.9*1.2,'ISO-E','FontSize',fontsize_legend);

h=legend('Constant','Time-varying','Location','best');
set(h,'fontsize',fontsize_legend);
set(gca,'linewidth',linewidth_gca,'fontsize',fontsize_gca)
set(gca,'GridLineStyle','--');
ylabel('PNCM - $$Q$$','Interpreter','latex','fontsize',fontsize_label)   % 纵坐标

% 设置输出保存图片的大小和格式
% hfig = figure(1);
% figWidth = 7.99;  % 设置图片宽度
% figHeight = 5;  % 设置图片高度
% set(hfig,'PaperUnits','centimeters'); % 图片尺寸所用单位
% set(hfig,'PaperPosition',[0 0 figWidth figHeight]);
% fileout = 'PNCM.'; % 输出图片的文件名
% print(hfig,[fileout,'tif'],'-r300','-dtiff'); % 设置图片格式、分辨率

%% PLOT(2)
figure(1)
hold on; box on; grid on;
set(gcf,'unit','centimeters','position',[10,10,8,5])    % 图形窗口在电脑屏幕上的位置和尺寸[左 下 宽 高]
% linewidth_line = 1.2;      % 图形线条宽度
linewidth_gca = 0.7;      % 横纵坐标轴宽度
fontsize_gca = 7;           % 横纵坐标轴刻度字体大小
fontsize_label = 9;         % 横纵坐标轴字体大小
fontsize_legend = 7;      % 图例字体大小


plot(R_cons(2,:),'LineWidth',2)
plot(R_true(2,:),'LineWidth',2)


line([501 501],[0 1.2],'linestyle','--','Color','k');
line([1001 1001],[0 1.2],'linestyle','--','Color','k');
text(200-10,0.9*1.2,'ISO-A','FontSize',fontsize_legend);
text(700-10,0.9*1.2,'ISO-C','FontSize',fontsize_legend);
text(1200-10,0.9*1.2,'ISO-E','FontSize',fontsize_legend);

h=legend('Constant','Time-varying','Location','best');
set(h,'fontsize',fontsize_legend);
set(gca,'linewidth',linewidth_gca,'fontsize',fontsize_gca)
set(gca,'GridLineStyle','--');
ylabel('PNCM - $$Q$$','Interpreter','latex','fontsize',fontsize_label)   % 纵坐标

% 设置输出保存图片的大小和格式
% hfig = figure(1);
% figWidth = 7.99;  % 设置图片宽度
% figHeight = 5;  % 设置图片高度
% set(hfig,'PaperUnits','centimeters'); % 图片尺寸所用单位
% set(hfig,'PaperPosition',[0 0 figWidth figHeight]);
% fileout = 'PNCM.'; % 输出图片的文件名
% print(hfig,[fileout,'tif'],'-r300','-dtiff'); % 设置图片格式、分辨率

%% PLOT(3)
figure(3)
hold on; box on; grid on;
set(gcf,'unit','centimeters','position',[10,10,8,5])    % 图形窗口在电脑屏幕上的位置和尺寸[左 下 宽 高]
% linewidth_line = 1.2;      % 图形线条宽度
linewidth_gca = 0.7;      % 横纵坐标轴宽度
fontsize_gca = 7;           % 横纵坐标轴刻度字体大小
fontsize_label = 9;         % 横纵坐标轴字体大小
fontsize_legend = 7;      % 图例字体大小

plot(R_cons(1,:),'LineWidth',2)
plot(R_cons(2,:),'LineWidth',2)
plot(R_cons(3,:),'LineWidth',2)
plot(R_true(1,:),'LineWidth',2)
plot(R_true(2,:),'LineWidth',2)
plot(R_true(3,:),'LineWidth',2)

line([501 501],[0 1.2],'linestyle','--','Color','k');
line([1001 1001],[0 1.2],'linestyle','--','Color','k');
text(200-10,0.9*1.2,'ISO-A','FontSize',fontsize_legend);
text(700-10,0.9*1.2,'ISO-C','FontSize',fontsize_legend);
text(1200-10,0.9*1.2,'ISO-E','FontSize',fontsize_legend);

h=legend('Constant','Time-varying','Location','best');
set(h,'fontsize',fontsize_legend);
set(gca,'linewidth',linewidth_gca,'fontsize',fontsize_gca)
set(gca,'GridLineStyle','--');
ylabel('PNCM - $$Q$$','Interpreter','latex','fontsize',fontsize_label)   % 纵坐标

% 设置输出保存图片的大小和格式
% hfig = figure(1);
% figWidth = 7.99;  % 设置图片宽度
% figHeight = 5;  % 设置图片高度
% set(hfig,'PaperUnits','centimeters'); % 图片尺寸所用单位
% set(hfig,'PaperPosition',[0 0 figWidth figHeight]);
% fileout = 'PNCM.'; % 输出图片的文件名
% print(hfig,[fileout,'tif'],'-r300','-dtiff'); % 设置图片格式、分辨率