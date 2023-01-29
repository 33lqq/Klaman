load fisheriris %下载样本数据集
indices = crossvalind('Kfold',species,10);%用了第一个句法，species为种类（个数）
species(1:50,1) = 1;
species(51:100,1) = 2;
species(101:150,1) = 3;
cp = classperf(species);%这里的<span style="font-family: Arial, Helvetica, sans-serif;">classperf是评估分类器性能的函数</span>
cp.ErrorRate
for i = 1:10
    test = (indices == i); train = ~test;%开始以第一个样本集为验证集，剩余的为训练集
%     class = classify(meas(test,:),meas(train,:),species(train,:));%用分类器测试
    T = ind2vec(species(train,:)');
    spread = 1.5;
    net = newpnn(meas(train,:)',T,spread);
    %使用测试集验证
%     testData = testData(:,setdiff(1:V, a)); %差分出非0索引
    y = net(meas(test,:)');
     ac = vec2ind(y);
    classperf(cp,net,test)%计算性能
    cp.ErrorRate
end
cp.ErrorRate
