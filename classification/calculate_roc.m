function [all_fpr,mean_tpr,Auc] = calculate_roc(xdata, ydata)
%UNTITLED7 此处提供此函数的摘要
%   此处提供详细说明
n_classes = 8;
% macro
% First aggregate all false positive rates
all_fpr = unique(cell2mat(xdata));

mean_tpr = zeros(size(all_fpr),'like',all_fpr);
for i=1:n_classes
    [xp,ia,~] = unique(xdata{i},'last','legacy') ;
    fp = ydata{i}(ia);
    mean_tpr = mean_tpr + interp1(xp, fp , all_fpr);
end
mean_tpr = mean_tpr / n_classes;

Auc = 0;
for i = 2:length(all_fpr)
    Auc = Auc + (mean_tpr(i) + mean_tpr(i-1))*(all_fpr(i) - all_fpr(i-1))/2;
end
end