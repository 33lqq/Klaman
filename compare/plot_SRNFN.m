%% RMSE
SRNFN_P1 = mean(SRNFN_P1, 3) .^ 0.25;        %% PECM
SRNFN_R1 = mean(SRNFN_R1, 3) .^ 0.25;        %% MNCM

SRNFN_P2 = mean(SRNFN_P2, 3) .^ 0.25; 
SRNFN_R2 = mean(SRNFN_R2, 3) .^ 0.25;

SRNFN_P3 = mean(SRNFN_P3, 3) .^ 0.25; 
SRNFN_R3 = mean(SRNFN_R3, 3) .^ 0.25;

SRNFN_P4 = mean(SRNFN_P4, 3) .^ 0.25;         %% PECM
SRNFN_P41 = mean(SRNFN_P41, 3) .^ 0.25;             %% PECM_with smooth
SRNFN_R4 = mean(SRNFN_R4, 3) .^ 0.25;         %% MNCM
SRNFN_R41 = mean(SRNFN_R41, 3) .^ 0.25;       %% MNCM_with smooth

%% ASRNFN
% ASRNFN_P1 = mean(mean(SRNFN_P1, 3)) .^ 0.25;  %% PECM
% ASRNFN_R1 = mean(mean(SRNFN_R1, 3)) .^ 0.25;  %% MNCM
% 
% ASRNFN_P2 = mean(mean(SRNFN_P2, 3)) .^ 0.25; 
% ASRNFN_R2 = mean(mean(SRNFN_R2, 3)) .^ 0.25;
% 
% ASRNFN_P3 = mean(mean(SRNFN_P3, 3)) .^ 0.25; 
% ASRNFN_R3 = mean(mean(SRNFN_R3, 3)) .^ 0.25;
% 
% ASRNFN_P4 = mean(mean(SRNFN_P4, 3)) .^ 0.25;   %% PECM
% ASRNFN_P41 = mean(mean(SRNFN_P41, 3)) .^ 0.25; %% PECM_with smooth
% ASRNFN_R4 = mean(mean(SRNFN_R4, 3)) .^ 0.25;   %% MNCM
% ASRNFN_R41 = mean(mean(SRNFN_R41, 3)) .^ 0.25; %% MNCM_with smooth
temp = mean(SRNFN_P1, 3);
Stretch1 = temp(2:501);
Stretch2 = temp(502:1001);
Stretch3 = temp(1002:1501);
ASRNFN_P1_1 = mean(Stretch1) .^ 0.25;   %% position
ASRNFN_P1_2 = mean(Stretch2) .^ 0.25;   %% position
ASRNFN_P1_3 = mean(Stretch3) .^ 0.25;   %% position
temp = mean(SRNFN_R1, 3);
Stretch1 = temp(2:501);
Stretch2 = temp(502:1001);
Stretch3 = temp(1002:1501);
ASRNFN_R1_1 = mean(Stretch1).^ 0.25; %% velocity
ASRNFN_R1_2 = mean(Stretch2).^ 0.25; %% velocity
ASRNFN_R1_3 = mean(Stretch3).^ 0.25; %% velocity

%
temp = mean(SRNFN_P2, 3);
Stretch1 = temp(2:501);
Stretch2 = temp(502:1001);
Stretch3 = temp(1002:1501);
ASRNFN_P2_1 = mean(Stretch1) .^ 0.25; 
ASRNFN_P2_2 = mean(Stretch2) .^ 0.25; 
ASRNFN_P2_3 = mean(Stretch3) .^ 0.25; 
temp = mean(SRNFN_R2, 3);
Stretch1 = temp(2:501);
Stretch2 = temp(502:1001);
Stretch3 = temp(1002:1501);
ASRNFN_R2_1 = mean(Stretch1) .^ 0.25;
ASRNFN_R2_2 = mean(Stretch2) .^ 0.25;
ASRNFN_R2_3 = mean(Stretch3) .^ 0.25;

%
temp = mean(SRNFN_P3, 3);
Stretch1 = temp(2:501);
Stretch2 = temp(502:1001);
Stretch3 = temp(1002:1501);
ASRNFN_P3_1 = mean(Stretch1) .^ 0.25; 
ASRNFN_P3_2 = mean(Stretch2) .^ 0.25; 
ASRNFN_P3_3 = mean(Stretch3) .^ 0.25; 
temp = mean(SRNFN_R3, 3);
Stretch1 = temp(2:501);
Stretch2 = temp(502:1001);
Stretch3 = temp(1002:1501);
ASRNFN_R3_1 = mean(Stretch1) .^ 0.25;
ASRNFN_R3_2 = mean(Stretch2) .^ 0.25;
ASRNFN_R3_3 = mean(Stretch3) .^ 0.25;

%
temp = mean(SRNFN_P41, 3);
Stretch1 = temp(2:501);
Stretch2 = temp(502:1001);
Stretch3 = temp(1002:1501);
ASRNFN_P41_1 = mean(Stretch1) .^ 0.25; 
ASRNFN_P41_2 = mean(Stretch2) .^ 0.25; 
ASRNFN_P41_3 = mean(Stretch3) .^ 0.25; 
temp = mean(SRNFN_R4, 3);
Stretch1 = temp(2:501);
Stretch2 = temp(502:1001);
Stretch3 = temp(1002:1501);
ASRNFN_R4_1 = mean(Stretch1) .^ 0.25;
ASRNFN_R4_2 = mean(Stretch2) .^ 0.25;
ASRNFN_R4_3 = mean(Stretch3) .^ 0.25;

%% PLOT(1)
subplot(1,2,1)
hold on; box on; grid on;
set(gcf,'unit','centimeters','position',[10,10,8,5])    % 图形窗口在电脑屏幕上的位置和尺寸[左 下 宽 高]
linewidth_line = 1.2;      % 图形线条宽度
linewidth_gca = 0.7;      % 横纵坐标轴宽度
fontsize_gca = 7;           % 横纵坐标轴刻度字体大小
fontsize_label = 9;         % 横纵坐标轴字体大小
fontsize_legend = 7;      % 图例字体大小
    
plot( SRNFN_P1,'linewidth',linewidth_line,'Color','#2878B5')
plot( SRNFN_P2,'linewidth',linewidth_line,'Color','#9AC9DB')
plot( SRNFN_P3,'linewidth',linewidth_line,'Color','#F8AC8C')
% plot( SRNFN_P4,'linewidth',linewidth_line)
plot( SRNFN_P41,'linewidth',linewidth_line,'Color','#C82423')
% h = legend('KFTCM','KFNCM','RSAKF','VBAKF-R','The proposed filter_withoutsmooth','The proposed filter');       % 图例


lim = 0.04;
line([501 501],[0 lim],'linestyle','--','Color','k');
line([1001 1001],[0 lim],'linestyle','--','Color','k');
text(200-20,0.9*lim,'ISO-A','FontSize',6.5);
text(700-20,0.9*lim,'ISO-C','FontSize',6.5);
text(1200-20,0.9*lim,'ISO-E','FontSize',6.5);

h = legend('KFNCM','RSAKF','VBAKF-R','The proposed filter','Location','west');       % 图例
set(h,'fontsize',fontsize_legend);
set(gca,'linewidth',linewidth_gca,'fontsize',fontsize_gca)
set(gca,'GridLineStyle','--');
ylabel('SRNFN of PECM','fontsize',fontsize_label)   % 纵坐标


%% PLOT(2)
subplot(1,2,2)
hold on; box on; grid on;
set(gcf,'unit','centimeters','position',[10,10,8,5])    % 图形窗口在电脑屏幕上的位置和尺寸[左 下 宽 高]
linewidth_line = 1.2;      % 图形线条宽度
linewidth_gca = 0.7;      % 横纵坐标轴宽度
fontsize_gca = 7;           % 横纵坐标轴刻度字体大小
fontsize_label = 9;         % 横纵坐标轴字体大小
fontsize_legend = 7;      % 图例字体大小

plot( SRNFN_R1,'linewidth',linewidth_line,'Color','#2878B5')
plot( SRNFN_R2,'linewidth',linewidth_line,'Color','#9AC9DB')
plot( SRNFN_R3,'linewidth',linewidth_line,'Color','#F8AC8C')
plot( SRNFN_R4,'linewidth',linewidth_line,'Color','#C82423')
% plot( SRNFN_R4    1,'Color','r','linewidth',linewidth_line)
% h = legend('KFTCM','KFNCM','RSAKF','VBAKF-R','The proposed filter_withoutsmooth','The proposed filter');       % 图例


lim =0.4;
line([501 501],[0 lim],'linestyle','--','Color','k');
line([1001 1001],[0 lim],'linestyle','--','Color','k');
text(200-10,0.9*lim,'ISO-A','FontSize',6.5);
text(700-10,0.9*lim,'ISO-C','FontSize',6.5);
text(1200-10,0.9*lim,'ISO-E','FontSize',6.5);

h = legend('KFNCM','RSAKF','VBAKF-R','The proposed filter','Location','west');       % 图例
set(h,'fontsize',fontsize_legend);
set(gca,'linewidth',linewidth_gca,'fontsize',fontsize_gca)
set(gca,'GridLineStyle','--');
ylabel('SRNFN of MNCM','fontsize',fontsize_label)   % 纵坐标

% 设置输出保存图片的大小和格式
hfig = figure(1);
figWidth = 16;  % 设置图片宽度
figHeight = 7.5;  % 设置图片高度
set(hfig,'PaperUnits','centimeters'); % 图片尺寸所用单位
set(hfig,'PaperPosition',[0 0 figWidth figHeight]);
fileout = 'SRNFN_11.'; % 输出图片的文件名
print(hfig,[fileout,'pdf'],'-r300','-dpdf'); % 设置图片格式、分辨率

