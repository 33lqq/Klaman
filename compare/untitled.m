
close all

Nsim=501;
for i=1
    figure(i)
    set(gcf,'unit','centimeters','position',[10,10,12,7.5])    % 图形窗口在电脑屏幕上的位置和尺寸[左 下 宽 高]
    linewidth_line = 1.5;      % 图形线条宽度
    markersize = 4;          % 图形标记点大小
    linewidth_gca = 0.7;      % 横纵坐标轴宽度
    fontsize_gca = 10;           % 横纵坐标轴刻度字体大小
    fontsize_label = 12;         % 横纵坐标轴字体大小
    fontsize_legend = 9;      % 图例字体大小

    plot(time(1:Nsim),R_cons(i,1:Nsim))
    hold on
%     plot(time(1:Nsim),y_lqr(i,1:Nsim))
    plot(time(1:Nsim),R_true(i,1:Nsim))
%     plot(time(1:Nsim),y_Empc(i,1:Nsim))

%     h=legend('Passive','LQR','MPC','EMPC');
%     h=legend('Passive','MPC','EMPC');
    h=legend('Constant','Time-varying');

    set(h,'fontsize',fontsize_legend);
    set(gca,'linewidth',linewidth_gca,'fontsize',fontsize_gca)
    set(gca,'GridLineStyle','--');
%     xlabel('Time-[s]','fontsize',fontsize_label) % 横坐标
%     ylabel({'Value'}, 'Interpreter','latex','fontsize',fontsize_label)   % 纵坐标

end
%% 保存图片
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
fileout = ['Time_varying','.']; % 输出图片的文件名
print(fig,[fileout,'pdf'],'-bestfit','-r300','-dpdf'); % 设置图片格式、分辨率


