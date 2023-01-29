function f = initialize_variables(N, M, V, min_range, max_range, data, label) %初始化种群
% f的大小→N×（V+M）前V个是维度，后两个是目标函数
min = min_range;
max = max_range; 
% digits(12)
K = M + V; %K是数组的总元素个数。为了便于计算，决策变量和目标函数串在一起形成一个数组。

% rand = unifrnd(0,1,N,V);
for i = 1 : N
    for j = 1 : V
        f(i,j) = min(j) + (max(j) - min(j))*rand(1); %f(i j)表示的是种群中第i个个体中的第j个决策变量，
        
        if f(i,j)>=0.5
            f(i,j) = 1;
        else
            f(i,j) = 0;
        end
%         if rand(i,j)>=0.5 
%             rand(i,j) = 1;
%     	else 
%             rand(i,j) = 0;
%         end
% 
%         f(i,j) =  rand(i,j);
        
        %这行代码为每个个体的所有决策变量在约束条件内随机取值,
    %怎么实现在区间内任意取值？
    end 
    f(i,V + 1: K) = evaluate_objective(f(i,:), M, V, data, label); % M是目标函数数量 V是决策变量个数
    %为了简化计算将对应的目标函数值储存在染色体的V + 1 到 M + V的位置。
end
end