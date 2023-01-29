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
%     view(-31,14) % 相机视角
    
    set(gca,'linewidth',linewidth_gca,'fontsize',fontsize_gca)
    set(gca,'GridLineStyle','--');
    set(gca,'dataaspectratio',[1 1 0.5],'projection','perspective','box','on')
    xlabel('Classification accuracy','fontsize',fontsize_label);
    ylabel('Inter-cluster distance','fontsize',fontsize_label);
    zlabel('Intra-cluster distance','fontsize',fontsize_label);
%     xlim([0.4 1])
    ylim([0.15 0.6])
%     zlim([a b])
    
    h = rotate3d;
    set(h, 'ActionPreCallback', 'set(gcf,''windowbuttonmotionfcn'',@align_axislabel)')
    set(h, 'ActionPostCallback', 'set(gcf,''windowbuttonmotionfcn'','''')')
    set(gcf, 'ResizeFcn', @align_axislabel)
    align_axislabel([], gca)
    axislabel_translation_slider;

    cb = colorbar;
    cb.Label.String = 'Number of feature subsets';
    
    hfig = figure(1);
%     figWidth = 12;  % 设置图片宽度
%     figHeight = 7.5;  % 设置图片高度
%     set(hfig,'PaperUnits','centimeters'); % 图片尺寸所用单位
%     set(hfig,'PaperPosition',[0 0 figWidth figHeight]);
   set(hfig,'PaperPositionMode','auto');
    fileout = 'naga2.'; % 输出图片的文件名
    % print(hfig,[fileout,'tif'],'-r300','-dtiff'); % 设置图片格式、分辨率
    print(hfig,[fileout,'pdf'],'-r600','-dpdf'); % 设置图片格式、分辨率
