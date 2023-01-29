%�Գ�ʼ��Ⱥ��ʼ���򣬿��ٷ�֧������
%ʹ�÷�֧���������Ⱥ�������򡣸ú�������ÿ�������Ӧ������ֵ��ӵ�����룬��һ�����еľ���
%��������ֵ��ӵ��������ӵ�Ⱦɫ�������
function f = non_domination_sort_mod(x, M, V) %��֧������
[N, ~] = size(x); %��N,V+M��,NΪ����x��������Ҳ����Ⱥ������
clear m
front = 1;
F(front).f = [];
individual = [];

for i = 1 : N
    individual(i).n = 0;%n�Ǹ���i��֧��ĸ�������
    individual(i).p = [];%p�Ǳ�����i֧��ĸ��弯��
    for j = 1 : N %��Ⱥ��
        dom_less = 0;
        dom_equal = 0;
        dom_more = 0;
        for k = 1 : M
            if (x(i,V + k) < x(j,V + k)) %�Ƚ���Ӧ�ȣ��Ƚ�������Ӧ��
                dom_less = dom_less + 1; %i���ж���jС��i֧��j
            elseif (x(i,V + k) == x(j,V + k)) % ����һ��
                dom_equal = dom_equal + 1;
            else
                dom_more = dom_more + 1;
            end
        end
        if dom_less == 0 && dom_equal ~= M %˵��i��j֧�䣬��Ӧ��n��1
            individual(i).n = individual(i).n + 1;
        elseif dom_more == 0 && dom_equal ~= M %˵��i֧��j����j����i��֧��ϼ���
            individual(i).p = [individual(i).p j];
        end
    end   
    if individual(i).n == 0 %����i��֧��ȼ�������ߣ�i����j��֧�䣩�����ڵ�ǰ���Ž⼯����Ӧ��Ⱦɫ����Я����������������Ϣ
        x(i,M + V + 1) = 1; %��ǰ���Ž�
        F(front).f = [F(front).f i]; %�ȼ�Ϊ1�ķ�֧�伯
    %���~=0����˵��i�ܱ�ĸ���֧�䣬���и���ȫ��ָ�궼����i
    end

end
%�ρ9�7�Ĵ�����Ϊ���ҳ��ȼ���9�0�ā9�6�9�6��⼯,�����һ�����ŵĸ��壬��33�͵���0������������i���Ͳ���Ҫ����33��
%�9�7�Ĵ�����Ϊ�˸�����������9�5�ּ�
while ~isempty(F(front).f)

   Q = [];%����9�2��front����
   for i = 1 : length(F(front).f)%ѭ����ǰ�9�6��⼯�еĸ���
       if ~isempty(individual(F(front).f(i)).p)%����i�Ё9�3�9�0���9�6��Ľ⼯
        	for j = 1 : length(individual(F(front).f(i)).p)%ѭ������i���9�6��⼯�еĸ���
            	individual(individual(F(front).f(i)).p(j)).n = ...%...��9�4�������9�2�9�5�����������ģ� ��9�7��9�4����j�ı��9�6�������1
                	individual(individual(F(front).f(i)).p(j)).n - 1;
        	   	if individual(individual(F(front).f(i)).p(j)).n == 0% ���q�ǁ9�6�9�6��⼯����Ł9�2����Q��
               		x(individual(F(front).f(i)).p(j),M + V + 1) = ...%����Ⱦ�9�0���мӁ9�2�ּ���Ϣ
                        front + 1;
                    Q = [Q individual(F(front).f(i)).p(j)];
                end
            end
       end
   end
   front =  front + 1;
   F(front).f = Q;
end

[temp,index_of_fronts] = sort(x(:,M + V + 1));%�Ը���Ĵ�������ȼ������������9�5�������� index_of_fronts��9�4������ֵ��Ӧԭ��������
for i = 1 : length(index_of_fronts)
    sorted_based_on_front(i,:) = x(index_of_fronts(i),:);%sorted_based_on_front�д�ŵ���x����������ȼ����������ľ���
end
current_index = 0;

%% Crowding distance ����ÿ�������ӵ����

for front = 1 : (length(F) - 1)%��9�7��1����Ϊ����55�9�5��9�7��F�����9�2��Ԫ��Ϊ�գ�������������ѭ�������ԁ9�2����length-1������ȼ�
    distance = 0;
    y = [];
    previous_index = current_index + 1;
    for i = 1 : length(F(front).f)
        y(i,:) = sorted_based_on_front(current_index + i,:);%y�д�ŵ�������ȼ�Ϊfront�ļ��Ͼ���
    end
    current_index = current_index + i;
    sorted_based_on_objective = [];%��Ż���ӵ����������ľ���
    for i = 1 : M
        [sorted_based_on_objective, index_of_objectives] = ...
            sort(y(:,V + i));%���Ձ9�0�꺯��ֵ����
        sorted_based_on_objective = [];
        for j = 1 : length(index_of_objectives)
            sorted_based_on_objective(j,:) = y(index_of_objectives(j),:);% sorted_based_on_objective��Ű��Ձ9�0�꺯��ֵ������x����
        end
        f_max = ...
            sorted_based_on_objective(length(index_of_objectives), V + i);%fmaxΪ�9�0�꺯����9�8ֵ fminΪ�9�0�꺯����9�3ֵ
        f_min = sorted_based_on_objective(1, V + i);
        y(index_of_objectives(length(index_of_objectives)),M + V + 1 + i)...%�������ĵځ9�2����������9�2������ľ�����Ϊ�9�2��9�8
            = Inf;
        y(index_of_objectives(1),M + V + 1 + i) = Inf;
         for j = 2 : length(index_of_objectives) - 1 %ѭ�������г��˵ځ9�2�������9�2���ĸ���,��ֻ��2�����壬�ⲿ�ֶ�û����
            next_obj  = sorted_based_on_objective(j + 1,V + i);
            previous_obj  = sorted_based_on_objective(j - 1,V + i);
            if (f_max - f_min == 0)
                y(index_of_objectives(j),M + V + 1 + i) = Inf;
            else
                y(index_of_objectives(j),M + V + 1 + i) = ...
                     (next_obj - previous_obj)/(f_max - f_min); %2050�ǵ�һ����Ӧ�ȵ�ӵ���ȣ�2051�ǵڶ�����Ӧ�ȵ�ӵ����
            end
         end
    end
    distance = [];
    distance(:,1) = zeros(length(F(front).f),1);
    for i = 1 : M
        distance(:,1) = distance(:,1) + y(:,M + V + 1 + i);
    end
    y(:,M + V + 2) = distance;
    y = y(:,1 : M + V + 2);
    z(previous_index:current_index,:) = y;
end
f = z(); %
