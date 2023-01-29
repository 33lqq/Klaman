%对初始种群开始排序，快速非支配排序
%使用非支配排序对种群进行排序。该函数返回每个个体对应的排序值和拥挤距离，是一个两列的矩阵
%并将排序值和拥挤距离添加到染色体矩阵中
function f = non_domination_sort_mod(x, M, V) %非支配排序
[N, ~] = size(x); %（N,V+M）,N为矩阵x的行数，也是种群的数量
clear m
front = 1;
F(front).f = [];
individual = [];

for i = 1 : N
    individual(i).n = 0;%n是个体i被支配的个体数量
    individual(i).p = [];%p是被个体i支配的个体集合
    for j = 1 : N %种群书
        dom_less = 0;
        dom_equal = 0;
        dom_more = 0;
        for k = 1 : M
            if (x(i,V + k) < x(j,V + k)) %比较适应度，比较两个适应度
                dom_less = dom_less + 1; %i所有都比j小，i支配j
            elseif (x(i,V + k) == x(j,V + k)) % 两个一样
                dom_equal = dom_equal + 1;
            else
                dom_more = dom_more + 1;
            end
        end
        if dom_less == 0 && dom_equal ~= M %说明i受j支配，相应的n加1
            individual(i).n = individual(i).n + 1;
        elseif dom_more == 0 && dom_equal ~= M %说明i支配j，把j加入i的支配合集中
            individual(i).p = [individual(i).p j];
        end
    end   
    if individual(i).n == 0 %个体i非支配等级排序最高（i不受j的支配），属于当前最优解集，相应的染色体中携带代表排序数的信息
        x(i,M + V + 1) = 1; %当前最优解
        F(front).f = [F(front).f i]; %等级为1的非支配集
    %如果~=0，则说明i受别的个体支配，另有个体全部指标都优于i
    end

end
%上⾯的代码是为了找出等级最⾼的⾮⽀配解集,如果有一个最优的个体，则33就等于0，对于其他的i，就不需要运行33了
%下⾯的代码是为了给其他个体进⾏分级
while ~isempty(F(front).f)

   Q = [];%存放下⼀个front集合
   for i = 1 : length(F(front).f)%循环当前⽀配解集中的个体
       if ~isempty(individual(F(front).f(i)).p)%个体i有⾃⼰所⽀配的解集
        	for j = 1 : length(individual(F(front).f(i)).p)%循环个体i所⽀配解集中的个体
            	individual(individual(F(front).f(i)).p(j)).n = ...%...表⽰的是与下⼀⾏代码是相连的， 这⾥表⽰个体j的被⽀配个数减1
                	individual(individual(F(front).f(i)).p(j)).n - 1;
        	   	if individual(individual(F(front).f(i)).p(j)).n == 0% 如果q是⾮⽀配解集，则放⼊集合Q中
               		x(individual(F(front).f(i)).p(j),M + V + 1) = ...%个体染⾊体中加⼊分级信息
                        front + 1;
                    Q = [Q individual(F(front).f(i)).p(j)];
                end
            end
       end
   end
   front =  front + 1;
   F(front).f = Q;
end

[temp,index_of_fronts] = sort(x(:,M + V + 1));%对个体的代表排序等级的列向量进⾏升序排序 index_of_fronts表⽰排序后的值对应原来的索引
for i = 1 : length(index_of_fronts)
    sorted_based_on_front(i,:) = x(index_of_fronts(i),:);%sorted_based_on_front中存放的是x矩阵按照排序等级升序排序后的矩阵
end
current_index = 0;

%% Crowding distance 计算每个个体的拥挤度

for front = 1 : (length(F) - 1)%这⾥减1是因为代码55⾏这⾥，F的最后⼀个元素为空，这样才能跳出循环。所以⼀共有length-1个排序等级
    distance = 0;
    y = [];
    previous_index = current_index + 1;
    for i = 1 : length(F(front).f)
        y(i,:) = sorted_based_on_front(current_index + i,:);%y中存放的是排序等级为front的集合矩阵
    end
    current_index = current_index + i;
    sorted_based_on_objective = [];%存放基于拥挤距离排序的矩阵
    for i = 1 : M
        [sorted_based_on_objective, index_of_objectives] = ...
            sort(y(:,V + i));%按照⽬标函数值排序
        sorted_based_on_objective = [];
        for j = 1 : length(index_of_objectives)
            sorted_based_on_objective(j,:) = y(index_of_objectives(j),:);% sorted_based_on_objective存放按照⽬标函数值排序后的x矩阵
        end
        f_max = ...
            sorted_based_on_objective(length(index_of_objectives), V + i);%fmax为⽬标函数最⼤值 fmin为⽬标函数最⼩值
        f_min = sorted_based_on_objective(1, V + i);
        y(index_of_objectives(length(index_of_objectives)),M + V + 1 + i)...%对排序后的第⼀个个体和最后⼀个个体的距离设为⽆穷⼤
            = Inf;
        y(index_of_objectives(1),M + V + 1 + i) = Inf;
         for j = 2 : length(index_of_objectives) - 1 %循环集合中除了第⼀个和最后⼀个的个体,当只有2个个体，这部分都没运行
            next_obj  = sorted_based_on_objective(j + 1,V + i);
            previous_obj  = sorted_based_on_objective(j - 1,V + i);
            if (f_max - f_min == 0)
                y(index_of_objectives(j),M + V + 1 + i) = Inf;
            else
                y(index_of_objectives(j),M + V + 1 + i) = ...
                     (next_obj - previous_obj)/(f_max - f_min); %2050是第一个适应度的拥挤度，2051是第二个适应度的拥挤度
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
