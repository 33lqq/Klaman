function f = evaluate_objective(x, M, V, data, label)

% f = [];
% f(1) = x(1);
% g = 1;
% sum = 0;
% for i = 2:V
%     sum = sum + x(i);
% end
% sum = 9*(sum / (V-1));
% g = g + sum;
% f(2) = g * (1 - sqrt(x(1) / g));
% end

f = [];
% M = 2; %目标函数数量，f1、f2
% V = 2146; %维度,Xi的变量个数，f(x) = x1 + x2 + ... + xV
%% 分类精度、特征数
%选取训练数据和测试数据
indices = crossvalind('Kfold',label,10); %数据集分成10份
for i=1:10
    % 测试集和测试标签
    test = (indices == i);
    testLabel = label(test , :); 
    testData = data(test , :);
    % 训练集和训练标签
    train = ~test;
    trainLabel = label(train , :);
    trainData = data(train , :);

    a = find(~x); %查找0索引
    trainData = trainData(:,setdiff(1:V, a)); %差分出非0索引
    %分类精度
    T = ind2vec(trainLabel');
    spread = 1.5;
    net = newpnn(trainData',T,spread);
    %使用测试集验证
    testData = testData(:,setdiff(1:V, a)); %差分出非0索引
    y = net(testData');
    ac = vec2ind(y);
    %统计分类准确率
    b = find(~(ac-testLabel')); %查找0值，即分类正确
    % length(b) %分类正确的数量
    Accuracy(i) = length(b)/length(testLabel'); %分类准确数
end

f(1) = 1 - mean(Accuracy);

%特征数
f(2) = 0;
for i = 1:V
    f(2) = f(2) + x(i);
end

% x= data' ,(2146×800)
%y = label' ,(1×800)
% Sw 类内散步矩阵，类内距离的平方形式
% Sb 类间散步矩阵
x = data(:,setdiff(1:V, a))'; % (2146×800)
y = label'; % (1×800)
[L,N] = size(x);  %有2146*800维
c = max(y);
%类内距离
% Sw
m=[];
Sw=zeros(1);
for i=1:1:c
    y_temp = (y==i);
    x_temp = x(:,y_temp);
    P(i)=sum(y_temp)/N;
    m(:,i)=(mean(x_temp'))';
    Sw=Sw+P(i)*cov(x_temp');  %矩阵形式
end
f(3) = trace(Sw);

%类间距离
% Sb
m0=(sum(((ones(L,1)*P).*m)'))';
Sb=zeros(1);
for i=1:c
    Sb=Sb+P(i)*((m(:,i)-m0)*(m(:,i)-m0)');  %矩阵形式
end
f(4) = -trace(Sb);





