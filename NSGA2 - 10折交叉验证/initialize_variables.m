function f = initialize_variables(N, M, V, min_range, max_range, data, label) %��ʼ����Ⱥ
% f�Ĵ�С��N����V+M��ǰV����ά�ȣ���������Ŀ�꺯��
min = min_range;
max = max_range; 
% digits(12)
K = M + V; %K���������Ԫ�ظ�����Ϊ�˱��ڼ��㣬���߱�����Ŀ�꺯������һ���γ�һ�����顣

% rand = unifrnd(0,1,N,V);
for i = 1 : N
    for j = 1 : V
        f(i,j) = min(j) + (max(j) - min(j))*rand(1); %f(i j)��ʾ������Ⱥ�е�i�������еĵ�j�����߱�����
        
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
        
        %���д���Ϊÿ����������о��߱�����Լ�����������ȡֵ,
    %��ôʵ��������������ȡֵ��
    end 
    f(i,V + 1: K) = evaluate_objective(f(i,:), M, V, data, label); % M��Ŀ�꺯������ V�Ǿ��߱�������
    %Ϊ�˼򻯼��㽫��Ӧ��Ŀ�꺯��ֵ������Ⱦɫ���V + 1 �� M + V��λ�á�
end
end