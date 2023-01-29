%% RMSE
RMSE = mean(RMS, 3) .^ 0.5;    %% position
RMSE1 = mean(RMS1, 3) .^ 0.5;  %% velocity

RMSE_1 = mean(RMS_1, 3) .^ 0.5; 
RMSE_1_1 = mean(RMS_1_1, 3) .^ 0.5;

RMSE_2 = mean(RMS_2, 3) .^ 0.5; 
RMSE_2_1 = mean(RMS_2_1, 3) .^ 0.5;

RMSE_3 = mean(RMS_3, 3) .^ 0.5; 
RMSE_3_1 = mean(RMS_3_1, 3) .^ 0.5;

RMSE_4 = mean(RMS_4, 3) .^ 0.5; 
RMSE_4_1 = mean(RMS_4_1, 3) .^ 0.5;

RMSE_41 = mean(RMS_41, 3) .^ 0.5; 
RMSE_41_1 = mean(RMS_41_1, 3) .^ 0.5;

%% ARMSE
% ARMSE = mean(mean(RMS, 3)) .^ 0.5;   %% position
% ARMSE1 = mean(mean(RMS1, 3) ).^ 0.5; %% velocity
% 
% ARMSE_1 = mean(mean(RMS_1, 3)) .^ 0.5; 
% ARMSE_1_1 = mean(mean(RMS_1_1, 3)) .^ 0.5;
% 
% ARMSE_2 = mean(mean(RMS_2, 3)) .^ 0.5; 
% ARMSE_2_1 = mean(mean(RMS_2_1, 3)) .^ 0.5;
% 
% ARMSE_3 = mean(mean(RMS_3, 3)) .^ 0.5; 
% ARMSE_3_1 = mean(mean(RMS_3_1, 3)) .^ 0.5;
% 
% ARMSE_4 = mean(mean(RMS_4, 3)) .^ 0.5; 
% ARMSE_4_1 = mean(mean(RMS_4_1, 3)) .^ 0.5;
% 
% ARMSE_41 = mean(mean(RMS_41, 3)) .^ 0.5; 
% ARMSE_41_1 = mean(mean(RMS_41_1, 3)) .^ 0.5;
% 
% ARMSE_51 = mean(RMSE_51.^ 2) .^ 0.5; 
% ARMSE_51_1 = mean(RMSE_51_1.^ 2) .^ 0.5;
%
temp = mean(RMS, 3);
Stretch1 = temp(2:501);
Stretch2 = temp(502:1001);
Stretch3 = temp(1002:1501);
ARMSE_1 = mean(Stretch1) .^ 0.5;   %% position
ARMSE_2 = mean(Stretch2) .^ 0.5;   %% position
ARMSE_3 = mean(Stretch3) .^ 0.5;   %% position
temp = mean(RMS1, 3);
Stretch1 = temp(2:501);
Stretch2 = temp(502:1001);
Stretch3 = temp(1002:1501);
ARMSE1_1 = mean(Stretch1).^ 0.5; %% velocity
ARMSE1_2 = mean(Stretch2).^ 0.5; %% velocity
ARMSE1_3 = mean(Stretch3).^ 0.5; %% velocity

%
temp = mean(RMS_1, 3);
Stretch1 = temp(2:501);
Stretch2 = temp(502:1001);
Stretch3 = temp(1002:1501);
ARMSE_1_1 = mean(Stretch1) .^ 0.5; 
ARMSE_1_2 = mean(Stretch2) .^ 0.5; 
ARMSE_1_3 = mean(Stretch3) .^ 0.5; 
temp = mean(RMS_1_1, 3);
Stretch1 = temp(2:501);
Stretch2 = temp(502:1001);
Stretch3 = temp(1002:1501);
ARMSE_1_1_1 = mean(Stretch1) .^ 0.5;
ARMSE_1_1_2 = mean(Stretch2) .^ 0.5;
ARMSE_1_1_3 = mean(Stretch3) .^ 0.5;

%
temp = mean(RMS_2, 3);
Stretch1 = temp(2:501);
Stretch2 = temp(502:1001);
Stretch3 = temp(1002:1501);
ARMSE_2_1 = mean(Stretch1) .^ 0.5; 
ARMSE_2_2 = mean(Stretch2) .^ 0.5; 
ARMSE_2_3 = mean(Stretch3) .^ 0.5; 
temp = mean(RMS_2_1, 3);
Stretch1 = temp(2:501);
Stretch2 = temp(502:1001);
Stretch3 = temp(1002:1501);
ARMSE_2_1_1 = mean(Stretch1) .^ 0.5;
ARMSE_2_1_2 = mean(Stretch2) .^ 0.5;
ARMSE_2_1_3 = mean(Stretch3) .^ 0.5;

%
temp = mean(RMS_3, 3);
Stretch1 = temp(2:501);
Stretch2 = temp(502:1001);
Stretch3 = temp(1002:1501);
ARMSE_3_1 = mean(Stretch1) .^ 0.5; 
ARMSE_3_2 = mean(Stretch2) .^ 0.5; 
ARMSE_3_3 = mean(Stretch3) .^ 0.5; 
temp = mean(RMS_3_1, 3);
Stretch1 = temp(2:501);
Stretch2 = temp(502:1001);
Stretch3 = temp(1002:1501);
ARMSE_3_1_1 = mean(Stretch1) .^ 0.5;
ARMSE_3_1_2 = mean(Stretch2) .^ 0.5;
ARMSE_3_1_3 = mean(Stretch3) .^ 0.5;

%
% temp = [RMSE_41(1:502),RMSE_4(503:end)].^2;
temp = mean(RMS_41, 3);
Stretch1 = temp(2:501);
Stretch2 = temp(502:1001);
Stretch3 = temp(1002:1501);
ARMSE_51_1 = mean(Stretch1) .^ 0.5; 
ARMSE_51_2 = mean(Stretch2) .^ 0.5; 
ARMSE_51_3 = mean(Stretch3) .^ 0.5;
% RMSE_51_1=[RMSE_41_1(1:501),RMSE_41_1(502:1001)+0.0055,RMSE_41_1(1002:1501)+0.0022];
% temp =RMSE_51_1.^ 2;
temp = mean(RMS_41_1, 3);
Stretch1 = temp(2:501);
Stretch2 = temp(502:1001);
Stretch3 = temp(1002:1501);
ARMSE_51_1_1 = mean(Stretch1) .^ 0.5;
ARMSE_51_1_2 = mean(Stretch2) .^ 0.5;
ARMSE_51_1_3 = mean(Stretch3) .^ 0.5;

%% PLOT(1)
subplot(1,2,1)
hold on; box on; grid on;
set(gcf,'unit','centimeters','position',[10,10,8,5])    % 图形窗口在电脑屏幕上的位置和尺寸[左 下 宽 高]
linewidth_line = 1.2;      % 图形线条宽度
linewidth_gca = 0.7;      % 横纵坐标轴宽度
fontsize_gca = 7;           % 横纵坐标轴刻度字体大小
fontsize_label = 9;         % 横纵坐标轴字体大小
fontsize_legend = 7;      % 图例字体大小

plot( RMSE,'linewidth',linewidth_line,'Color','#2878B5')
plot( RMSE_1,'linewidth',linewidth_line,'Color','#9AC9DB')
plot( RMSE_2,'linewidth',linewidth_line,'Color','#F8AC8C')
plot( RMSE_3,'linewidth',linewidth_line,'Color','#FF8884')
plot( RMSE_41,'linewidth',linewidth_line,'Color','#C82423')
% RMSE_5=[RMSE_41(1:502),RMSE_4(503:end)];
% plot( RMSE_5,'linewidth',linewidth_line,'Color','#C82423')
% plot( RMSE_41,'linewidth',linewidth_line,'Color','#C82423')
% h = legend('KFTCM','KFNCM','RSAKF','VBAKF-R','The proposed filter_withoutsmooth','The proposed filter');       % 图例


lim=6e-3;
line([501 501],[0 lim],'linestyle','--','Color','k');
line([1001 1001],[0 lim],'linestyle','--','Color','k');
text(200-20,0.9*lim,'ISO-E','FontSize',6.5);
text(700-20,0.9*lim,'ISO-C','FontSize',6.5);
text(1200-20,0.9*lim,'ISO-A','FontSize',6.5);

h = legend('KFTCM','KFNCM','RSAKF','VBAKF-R','The proposed filter','Location','east');       % 图例
set(h,'fontsize',fontsize_legend);
set(gca,'linewidth',linewidth_gca,'fontsize',fontsize_gca)
set(gca,'GridLineStyle','--');
ylabel('RMSE_{pos}-[m]','fontsize',fontsize_label)   % 纵坐标


%% PLOT(2)
subplot(1,2,2)
hold on; box on; grid on;
set(gcf,'unit','centimeters','position',[10,10,8,5])    % 图形窗口在电脑屏幕上的位置和尺寸[左 下 宽 高]
linewidth_line = 1.2;      % 图形线条宽度
linewidth_gca = 0.7;      % 横纵坐标轴宽度
fontsize_gca = 7;           % 横纵坐标轴刻度字体大小
fontsize_label = 9;         % 横纵坐标轴字体大小
fontsize_legend = 7;      % 图例字体大小

plot( RMSE1,'linewidth',linewidth_line,'Color','#2878B5')
plot( RMSE_1_1,'linewidth',linewidth_line,'Color','#9AC9DB')
plot( RMSE_2_1,'linewidth',linewidth_line,'Color','#F8AC8C')
plot( RMSE_3_1,'linewidth',linewidth_line,'Color','#FF8884')
% plot( RMSE_4_1,'linewidth',linewidth_line,'Color','#C82423')
plot( RMSE_41_1,'linewidth',linewidth_line,'Color','#C82423')
% RMSE_51_1=[RMSE_41_1(1:501),RMSE_41_1(502:1001)+0.0055,RMSE_41_1(1002:1501)+0.0022];
% plot( RMSE_51_1,'linewidth',linewidth_line,'Color','#C82423')
% h = legend('KFTCM','KFNCM','RSAKF','VBAKF-R','The proposed filter_withoutsmooth','The proposed filter');       % 图例


lim =0.08;
line([501 501],[0 lim],'linestyle','--','Color','k');
line([1001 1001],[0 lim],'linestyle','--','Color','k');
text(200-20,0.9*lim,'ISO-E','FontSize',6.5);
text(700-20,0.9*lim,'ISO-C','FontSize',6.5);
text(1200-20,0.9*lim,'ISO-A','FontSize',6.5);

h = legend('KFTCM','KFNCM','RSAKF','VBAKF-R','The proposed filter','Location','northeast');       % 图例
set(h,'fontsize',fontsize_legend);
set(gca,'linewidth',linewidth_gca,'fontsize',fontsize_gca)
set(gca,'GridLineStyle','--');
ylabel('RMSE_{vel}-[m/s]','fontsize',fontsize_label)   % 纵坐标

% 设置输出保存图片的大小和格式
hfig = figure(1);
figWidth = 16;  % 设置图片宽度
figHeight = 7.5;  % 设置图片高度 %10
set(hfig,'PaperUnits','centimeters'); % 图片尺寸所用单位
set(hfig,'PaperPosition',[0 0 figWidth figHeight]);
fileout = 'RMSE_without_1.'; % 输出图片的文件名
print(hfig,[fileout,'pdf'],'-r300','-dpdf'); % 设置图片格式、分辨率



