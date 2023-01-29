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
% M = 2; %Ŀ�꺯��������f1��f2
% V = 2146; %ά��,Xi�ı���������f(x) = x1 + x2 + ... + xV
%% ���ྫ�ȡ�������
%ѡȡѵ�����ݺͲ�������
indices = crossvalind('Kfold',label,10); %���ݼ��ֳ�10��
for i=1:10
    % ���Լ��Ͳ��Ա�ǩ
    test = (indices == i);
    testLabel = label(test , :); 
    testData = data(test , :);
    % ѵ������ѵ����ǩ
    train = ~test;
    trainLabel = label(train , :);
    trainData = data(train , :);

    a = find(~x); %����0����
    trainData = trainData(:,setdiff(1:V, a)); %��ֳ���0����
    %���ྫ��
    T = ind2vec(trainLabel');
    spread = 1.5;
    net = newpnn(trainData',T,spread);
    %ʹ�ò��Լ���֤
    testData = testData(:,setdiff(1:V, a)); %��ֳ���0����
    y = net(testData');
    ac = vec2ind(y);
    %ͳ�Ʒ���׼ȷ��
    b = find(~(ac-testLabel')); %����0ֵ����������ȷ
    % length(b) %������ȷ������
    Accuracy(i) = length(b)/length(testLabel'); %����׼ȷ��
end

f(1) = 1 - mean(Accuracy);

%������
f(2) = 0;
for i = 1:V
    f(2) = f(2) + x(i);
end

% x= data' ,(2146��800)
%y = label' ,(1��800)
% Sw ����ɢ���������ھ����ƽ����ʽ
% Sb ���ɢ������
x = data(:,setdiff(1:V, a))'; % (2146��800)
y = label'; % (1��800)
[L,N] = size(x);  %��2146*800ά
c = max(y);
%���ھ���
% Sw
m=[];
Sw=zeros(1);
for i=1:1:c
    y_temp = (y==i);
    x_temp = x(:,y_temp);
    P(i)=sum(y_temp)/N;
    m(:,i)=(mean(x_temp'))';
    Sw=Sw+P(i)*cov(x_temp');  %������ʽ
end
f(3) = trace(Sw);

%������
% Sb
m0=(sum(((ones(L,1)*P).*m)'))';
Sb=zeros(1);
for i=1:c
    Sb=Sb+P(i)*((m(:,i)-m0)*(m(:,i)-m0)');  %������ʽ
end
f(4) = -trace(Sb);





