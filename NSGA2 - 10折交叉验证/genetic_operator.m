%交叉算法选择的是单点交叉，变异算法选择的是多项式变异
function f  = genetic_operator(parent_chromosome, M, V, mu, mum, l_limit, u_limit, data, label) % l是最小值，u是最大值
[N,m] = size(parent_chromosome); %N是交配池中的个体数量，种群数

clear m
p = 1;
was_crossover = 0;%是否交叉标志位
was_mutation = 0;%是否变异标志位


for i = 1 : N %这⾥虽然循环N次，但是每次循环都会有概率产⽣2个或者1个⼦代，所以最终产⽣的⼦代个体数量⼤约是2N个
    if rand(1) < 0.9%交叉概率0.9,高概率
        child_1 = [];
        child_2 = [];
        parent_1 = round(N*rand(1));
        if parent_1 < 1
            parent_1 = 1;
        end
        parent_2 = round(N*rand(1));
        if parent_2 < 1
            parent_2 = 1;
        end
        while isequal(parent_chromosome(parent_1,:),parent_chromosome(parent_2,:))
            parent_2 = round(N*rand(1));
            if parent_2 < 1
                parent_2 = 1;
            end
        end
        %用于交叉操作的两个父代种群
        parent_1 = parent_chromosome(parent_1,:);
        parent_2 = parent_chromosome(parent_2,:);

        %生成子代元素，单点交叉
        j = randi([1,V-1]);
        child_1 = ...
            [parent_1(1:j) parent_2(j+1:V)];
        child_2 = ...
             [parent_2(1:j) parent_1(j+1:V)];

        child_1(:,V + 1: M + V) = evaluate_objective(child_1, M, V, data, label);
        child_2(:,V + 1: M + V) = evaluate_objective(child_2, M, V, data, label);
        was_crossover = 1;
        was_mutation = 0;
    else %变异操作,随机多点位翻转突变
        parent_3 = round(N*rand(1));
        if parent_3 < 1
            parent_3 = 1;
        end
        child_3 = parent_chromosome(parent_3,:);
        %随机n个点变异操作
        n = randi([1,V]);
        s = randi(V,1,n);
        for j = 1:n
            if child_3(s(j)) == 0
                child_3(s(j)) = 1;
            elseif child_3(s(j)) == 1
                child_3(s(j)) = 0;
            end
        end

        child_3(:,V + 1: M + V) = evaluate_objective(child_3, M, V, data, label);
        was_mutation = 1;
        was_crossover = 0;
    end
    if was_crossover
        child(p,:) = child_1;
        child(p+1,:) = child_2;
        was_cossover = 0;
        p = p + 2;
    elseif was_mutation
        child(p,:) = child_3(1,1 : M + V);
        was_mutation = 0;
        p = p + 1;
    end
end
f = child;
end
