%% naive bayes
openfig('bayes_1.fig');  % 读取图片1
obj=get(gca,'Children');  % 固定读取方法，我不知道原理，也无需深究，读取全部                       % 然后选中最后一个变量即为x，y值
len_obj=length(obj);      % 获取长度，方便访问最后一个变量
xdata_bayes{:,1}=get(obj(len_obj),'Xdata');   %分别获取x，y
ydata_bayes{:,1}=get(obj(len_obj),'Ydata');

openfig('bayes_2.fig')    %  Figure2与Figure1同理
obj2=get(gca,'children');
len_obj2=length(obj2);
xdata_bayes{:,2}=get(obj2(len_obj2),'xdata');
ydata_bayes{:,2}=get(obj2(len_obj2),'ydata');

openfig('bayes_3.fig');
obj3=get(gca,'children');
len_obj3=length(obj3);
xdata_bayes{:,3}=get(obj3(len_obj3),'xdata');
ydata_bayes{:,3}=get(obj3(len_obj3),'ydata');

openfig('bayes_4.fig');
obj4=get(gca,'children');
len_obj4=length(obj4);
xdata_bayes{:,4}=get(obj4(len_obj4),'xdata');
ydata_bayes{:,4}=get(obj4(len_obj4),'ydata');

openfig('bayes_5.fig');
obj5=get(gca,'children');
len_obj5=length(obj5);
xdata_bayes{:,5}=get(obj5(len_obj5),'xdata');
ydata_bayes{:,5}=get(obj5(len_obj5),'ydata');

openfig('bayes_6.fig');
obj6=get(gca,'children');
len_obj6=length(obj6);
xdata_bayes{:,6}=get(obj6(len_obj6),'xdata');
ydata_bayes{:,6}=get(obj6(len_obj6),'ydata');

openfig('bayes_7.fig');
obj7=get(gca,'children');
len_obj7=length(obj7);
xdata_bayes{:,7}=get(obj7(len_obj7),'xdata');
ydata_bayes{:,7}=get(obj7(len_obj7),'ydata');

openfig('bayes_8.fig');
obj8=get(gca,'children');
len_obj8=length(obj8);
xdata_bayes{:,8}=get(obj8(len_obj8),'xdata');
ydata_bayes{:,8}=get(obj8(len_obj8),'ydata');

close("all")             % 关闭全部窗口非常重要
clearvars -except ydata_bayes xdata_bayes

%% tree
openfig('tree_1.fig');  % 读取图片1
obj=get(gca,'Children');  % 固定读取方法，我不知道原理，也无需深究，读取全部                       % 然后选中最后一个变量即为x，y值
len_obj=length(obj);      % 获取长度，方便访问最后一个变量
xdata_tree{:,1}=get(obj(len_obj),'Xdata');   %分别获取x，y
ydata_tree{:,1}=get(obj(len_obj),'Ydata');

openfig('tree_2.fig')    %  Figure2与Figure1同理
obj2=get(gca,'children');
len_obj2=length(obj2);
xdata_tree{:,2}=get(obj2(len_obj2),'xdata');
ydata_tree{:,2}=get(obj2(len_obj2),'ydata');

openfig('tree_3.fig');
obj3=get(gca,'children');
len_obj3=length(obj3);
xdata_tree{:,3}=get(obj3(len_obj3),'xdata');
ydata_tree{:,3}=get(obj3(len_obj3),'ydata');

openfig('tree_4.fig');
obj4=get(gca,'children');
len_obj4=length(obj4);
xdata_tree{:,4}=get(obj4(len_obj4),'xdata');
ydata_tree{:,4}=get(obj4(len_obj4),'ydata');

openfig('tree_5.fig');
obj5=get(gca,'children');
len_obj5=length(obj5);
xdata_tree{:,5}=get(obj5(len_obj5),'xdata');
ydata_tree{:,5}=get(obj5(len_obj5),'ydata');

openfig('tree_6.fig');
obj6=get(gca,'children');
len_obj6=length(obj6);
xdata_tree{:,6}=get(obj6(len_obj6),'xdata');
ydata_tree{:,6}=get(obj6(len_obj6),'ydata');

openfig('tree_7.fig');
obj7=get(gca,'children');
len_obj7=length(obj7);
xdata_tree{:,7}=get(obj7(len_obj7),'xdata');
ydata_tree{:,7}=get(obj7(len_obj7),'ydata');

openfig('tree_8.fig');
obj8=get(gca,'children');
len_obj8=length(obj8);
xdata_tree{:,8}=get(obj8(len_obj8),'xdata');
ydata_tree{:,8}=get(obj8(len_obj8),'ydata');

close("all")             % 关闭全部窗口非常重要
clearvars -except ydata_bayes xdata_bayes ...
    xdata_tree ydata_tree
%% cnn
openfig('cnn_1.fig');  % 读取图片1
obj=get(gca,'Children');  % 固定读取方法，我不知道原理，也无需深究，读取全部                       % 然后选中最后一个变量即为x，y值
len_obj=length(obj);      % 获取长度，方便访问最后一个变量
xdata_cnn{:,1}=get(obj(len_obj),'Xdata');   %分别获取x，y
ydata_cnn{:,1}=get(obj(len_obj),'Ydata');

openfig('cnn_2.fig')    %  Figure2与Figure1同理
obj2=get(gca,'children');
len_obj2=length(obj2);
xdata_cnn{:,2}=get(obj2(len_obj2),'xdata');
ydata_cnn{:,2}=get(obj2(len_obj2),'ydata');

openfig('cnn_3.fig');
obj3=get(gca,'children');
len_obj3=length(obj3);
xdata_cnn{:,3}=get(obj3(len_obj3),'xdata');
ydata_cnn{:,3}=get(obj3(len_obj3),'ydata');

openfig('cnn_4.fig');
obj4=get(gca,'children');
len_obj4=length(obj4);
xdata_cnn{:,4}=get(obj4(len_obj4),'xdata');
ydata_cnn{:,4}=get(obj4(len_obj4),'ydata');

openfig('cnn_5.fig');
obj5=get(gca,'children');
len_obj5=length(obj5);
xdata_cnn{:,5}=get(obj5(len_obj5),'xdata');
ydata_cnn{:,5}=get(obj5(len_obj5),'ydata');

openfig('cnn_6.fig');
obj6=get(gca,'children');
len_obj6=length(obj6);
xdata_cnn{:,6}=get(obj6(len_obj6),'xdata');
ydata_cnn{:,6}=get(obj6(len_obj6),'ydata');

openfig('cnn_7.fig');
obj7=get(gca,'children');
len_obj7=length(obj7);
xdata_cnn{:,7}=get(obj7(len_obj7),'xdata');
ydata_cnn{:,7}=get(obj7(len_obj7),'ydata');

openfig('cnn_8.fig');
obj8=get(gca,'children');
len_obj8=length(obj8);
xdata_cnn{:,8}=get(obj8(len_obj8),'xdata');
ydata_cnn{:,8}=get(obj8(len_obj8),'ydata');

close("all")             % 关闭全部窗口非常重要
clearvars -except ydata_bayes xdata_bayes ...
    xdata_tree ydata_tree...
    xdata_cnn ydata_cnn

%% knn
openfig('knn_1.fig');  % 读取图片1
obj=get(gca,'Children');  % 固定读取方法，我不知道原理，也无需深究，读取全部                       % 然后选中最后一个变量即为x，y值
len_obj=length(obj);      % 获取长度，方便访问最后一个变量
xdata_knn{:,1}=get(obj(len_obj),'Xdata');   %分别获取x，y
ydata_knn{:,1}=get(obj(len_obj),'Ydata');

openfig('knn_2.fig')    %  Figure2与Figure1同理
obj2=get(gca,'children');
len_obj2=length(obj2);
xdata_knn{:,2}=get(obj2(len_obj2),'xdata');
ydata_knn{:,2}=get(obj2(len_obj2),'ydata');

openfig('knn_3.fig');
obj3=get(gca,'children');
len_obj3=length(obj3);
xdata_knn{:,3}=get(obj3(len_obj3),'xdata');
ydata_knn{:,3}=get(obj3(len_obj3),'ydata');

openfig('knn_4.fig');
obj4=get(gca,'children');
len_obj4=length(obj4);
xdata_knn{:,4}=get(obj4(len_obj4),'xdata');
ydata_knn{:,4}=get(obj4(len_obj4),'ydata');

openfig('knn_5.fig');
obj5=get(gca,'children');
len_obj5=length(obj5);
xdata_knn{:,5}=get(obj5(len_obj5),'xdata');
ydata_knn{:,5}=get(obj5(len_obj5),'ydata');

openfig('knn_6.fig');
obj6=get(gca,'children');
len_obj6=length(obj6);
xdata_knn{:,6}=get(obj6(len_obj6),'xdata');
ydata_knn{:,6}=get(obj6(len_obj6),'ydata');

openfig('knn_7.fig');
obj7=get(gca,'children');
len_obj7=length(obj7);
xdata_knn{:,7}=get(obj7(len_obj7),'xdata');
ydata_knn{:,7}=get(obj7(len_obj7),'ydata');

openfig('knn_8.fig');
obj8=get(gca,'children');
len_obj8=length(obj8);
xdata_knn{:,8}=get(obj8(len_obj8),'xdata');
ydata_knn{:,8}=get(obj8(len_obj8),'ydata');

close("all")             % 关闭全部窗口非常重要
clearvars -except ydata_bayes xdata_bayes ...
    xdata_tree ydata_tree...
    xdata_cnn ydata_cnn...
    xdata_knn ydata_knn

%% svm
openfig('svm_1.fig');  % 读取图片1
obj=get(gca,'Children');  % 固定读取方法，我不知道原理，也无需深究，读取全部                       % 然后选中最后一个变量即为x，y值
len_obj=length(obj);      % 获取长度，方便访问最后一个变量
xdata_svm{:,1}=get(obj(len_obj),'Xdata');   %分别获取x，y
ydata_svm{:,1}=get(obj(len_obj),'Ydata');

openfig('svm_2.fig')    %  Figure2与Figure1同理
obj2=get(gca,'children');
len_obj2=length(obj2);
xdata_svm{:,2}=get(obj2(len_obj2),'xdata');
ydata_svm{:,2}=get(obj2(len_obj2),'ydata');

openfig('svm_3.fig');
obj3=get(gca,'children');
len_obj3=length(obj3);
xdata_svm{:,3}=get(obj3(len_obj3),'xdata');
ydata_svm{:,3}=get(obj3(len_obj3),'ydata');

openfig('svm_4.fig');
obj4=get(gca,'children');
len_obj4=length(obj4);
xdata_svm{:,4}=get(obj4(len_obj4),'xdata');
ydata_svm{:,4}=get(obj4(len_obj4),'ydata');

openfig('svm_5.fig');
obj5=get(gca,'children');
len_obj5=length(obj5);
xdata_svm{:,5}=get(obj5(len_obj5),'xdata');
ydata_svm{:,5}=get(obj5(len_obj5),'ydata');

openfig('svm_6.fig');
obj6=get(gca,'children');
len_obj6=length(obj6);
xdata_svm{:,6}=get(obj6(len_obj6),'xdata');
ydata_svm{:,6}=get(obj6(len_obj6),'ydata');

openfig('svm_7.fig');
obj7=get(gca,'children');
len_obj7=length(obj7);
xdata_svm{:,7}=get(obj7(len_obj7),'xdata');
ydata_svm{:,7}=get(obj7(len_obj7),'ydata');

openfig('svm_8.fig');
obj8=get(gca,'children');
len_obj8=length(obj8);
xdata_svm{:,8}=get(obj8(len_obj8),'xdata');
ydata_svm{:,8}=get(obj8(len_obj8),'ydata');

close("all")             % 关闭全部窗口非常重要
clearvars -except ydata_bayes xdata_bayes ...
    xdata_tree ydata_tree...
    xdata_cnn ydata_cnn...
    xdata_knn ydata_knn...
    xdata_svm ydata_svm

%% plot
% 接下来只需要正常画图即可 
figure(1)
set(gcf,'unit','centimeters','position',[10,10,8,5])    % 图形窗口在电脑屏幕上的位置和尺寸[左 下 宽 高]
linewidth_line = 1.2;      % 图形线条宽度
% markersize = 2.5;          % 图形标记点大小
linewidth_gca = 0.7;      % 横纵坐标轴宽度
fontsize_gca = 7;           % 横纵坐标轴刻度字体大小
fontsize_label = 9;         % 横纵坐标轴字体大小
fontsize_legend = 7;      % 图例字体大小

[fpr_bayes, tpr_bayes, Auc_bayes] = calculate_roc(xdata_bayes, ydata_bayes);
plot(fpr_bayes, tpr_bayes,'linewidth',linewidth_line)
hold on;  
[fpr_cnn, tpr_cnn, Auc_cnn] = calculate_roc(xdata_cnn, ydata_cnn);
plot(fpr_cnn, tpr_cnn,'linewidth',linewidth_line)

[fpr_knn, tpr_knn, Auc_knn] = calculate_roc(xdata_knn, ydata_knn);
plot(fpr_knn, tpr_knn,'linewidth',linewidth_line)

[fpr_svm, tpr_svm, Auc_svm] = calculate_roc(xdata_svm, ydata_svm);
plot(fpr_svm, tpr_svm,'linewidth',linewidth_line)

[fpr_tree, tpr_tree, Auc_tree] = calculate_roc(xdata_tree, ydata_tree);
plot(fpr_tree, tpr_tree,'linewidth',linewidth_line)

xlim([0 1])           % X轴坐标范围
ylim([0 1.05])       % Y轴坐标范围
% ['NB (area=:'
h = legend(['NB (area = ',num2str(Auc_bayes,'%0.4f'),')'],...
    ['N-CNN (area = ',num2str(Auc_cnn,'%0.4f'),')'],...
    ['KNN (area = ',num2str(Auc_knn,'%0.4f'),')'],...
    ['L-SVM (area = ',num2str(Auc_svm,'%0.4f'),')'],...
    ['DT (area = ',num2str(Auc_tree,'%0.4f'),')'],'Location','southeast');       % 图例
set(h,'fontsize',fontsize_legend);
set(gca,'linewidth',linewidth_gca,'fontsize',fontsize_gca)
set(gca,'GridLineStyle','--');
xlabel('False Positive Rate','fontsize',fontsize_label) % 横坐标
ylabel('True Positive Rate','fontsize',fontsize_label)   % 纵坐标

% 设置输出保存图片的大小和格式
hfig = figure(1);
figWidth = 7.99;  % 设置图片宽度
figHeight = 5;  % 设置图片高度
set(hfig,'PaperUnits','centimeters'); % 图片尺寸所用单位
set(hfig,'PaperPosition',[0 0 figWidth figHeight]);
fileout = 'multi-calss ROC.'; % 输出图片的文件名
% print(hfig,[fileout,'tif'],'-r300','-dtiff'); % 设置图片格式、分辨率
print(hfig,[fileout,'eps'],'-r300','-deps'); % 设置图片格式、分辨率