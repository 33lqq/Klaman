%�����㷨ѡ����ǵ��㽻�棬�����㷨ѡ����Ƕ���ʽ����
function f  = genetic_operator(parent_chromosome, M, V, mu, mum, l_limit, u_limit, data, label) % l����Сֵ��u�����ֵ
[N,m] = size(parent_chromosome); %N�ǽ�����еĸ�����������Ⱥ��

clear m
p = 1;
was_crossover = 0;%�Ƿ񽻲��־λ
was_mutation = 0;%�Ƿ�����־λ


for i = 1 : N %��9�7��Ȼѭ��N�Σ�����ÿ��ѭ�������и��ʲ��9�12������1���9�0�����������ղ��9�1�ā9�0�����������9�8Լ��2N��
    if rand(1) < 0.9%�������0.9,�߸���
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
        %���ڽ������������������Ⱥ
        parent_1 = parent_chromosome(parent_1,:);
        parent_2 = parent_chromosome(parent_2,:);

        %�����Ӵ�Ԫ�أ����㽻��
        j = randi([1,V-1]);
        child_1 = ...
            [parent_1(1:j) parent_2(j+1:V)];
        child_2 = ...
             [parent_2(1:j) parent_1(j+1:V)];

        child_1(:,V + 1: M + V) = evaluate_objective(child_1, M, V, data, label);
        child_2(:,V + 1: M + V) = evaluate_objective(child_2, M, V, data, label);
        was_crossover = 1;
        was_mutation = 0;
    else %�������,������λ��תͻ��
        parent_3 = round(N*rand(1));
        if parent_3 < 1
            parent_3 = 1;
        end
        child_3 = parent_chromosome(parent_3,:);
        %���n����������
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
