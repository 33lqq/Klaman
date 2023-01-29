%****************************************
%          特征提取函数
%****************************************
function result = main_feature()
sample_num = 10000*15; %每一个文件是1500长，每10个点作为1个data
sensor_num = 3; %参数变量，y1-y3
feature_num = (11+64)*3; %特征量，3×12,1是速度
% num = 65;

for number_index =1:3
    number = number_index;
    level = ['A','B','C','D','E','F','G','H'];

    disp('-------------------------')
    disp('       正在提取数据       ')
    disp('-------------------------')

    data = zeros(sample_num, feature_num, 'double');
    %     fprintf('Processing E set: %03d\n', sample_index )
    %     disp(['Processing',' ',level(number), ' ','set:', num2str(sample_index )])
    %遍历文件夹文件
    maindir = 'D:/MATLAB/kalman/feature';
    dirout = dir(fullfile(maindir));
    % 剔除 ./..
    for k = length(dirout):-1:1
        if sum(strcmp(dirout(k).name,{'.','..'}))
            dirout(k) = [];
        end
    end
    % 筛选指定
    str = level(number);
    for k = length(dirout):-1:1
        if ~contains(dirout(k).name,str)
            dirout(k) = [];
        end
    end
    csv = [];
    for k = 1:length(dirout)
        subdir = ['D:/MATLAB/kalman/feature/',dirout(k).name];
        read_path = fullfile(subdir,'*.csv'); %合并a，b两个文件名,会重复
        csv = cat(1,csv,dir(read_path)) ;%在这个子文件夹下找后缀为csv的文件
    end
    j=1;
    for sample_index = 1:10000
        csv(sample_index).name;

        a = strsplit(csv(sample_index).name,{'.c'}); %切分得到原始时域名
        a = strfind({csv.name},a(1)); %1的位置就是

        b = 1;
        f = @(x) isequaln(x,b);
        inedx = find(cellfun(f,a));


        for i = 1:length(inedx)
            csvpath = fullfile(csv(inedx(i)).folder,csv(inedx(i)).name);
            if i == 1
                pre_data = readmatrix(csvpath); %1500×3

                for num = 1:15       % 从[1-10],[11-20],...总共15个数据
                    data_temp = pre_data(1+(num-1)*100:num*100,:);
                    k=0;
                    for m = 1:sensor_num %3
                        k=k+1;
                        data(num+(j-1)*15  , 1+(k-1)*(11+64):k*(11+64)) = extract_fea_ori(data_temp(:,m));
                    end

                end
                j=j+1;

            end
        end

    end
    result = data;

    currentRoadLevel = sprintf(['D:/MATLAB/kalman/feature/','%s.mat'],level(number));
    save(currentRoadLevel,'result')
    disp('-------------------------')
    disp('       完成提取，已自动保存     ')
    disp('-------------------------')
end

%% 时域
function fea_time = time_fea(a) %均方根
fea_1 = rms(a); %均方根
fea_2 = var(a); %方差
fea_3 = max(a); %最大值
fea_4 = skewness(a); %偏度
fea_5 = kurtosis(a); %峰度
fea_6 = max(a) - min(a); %峰峰值
fea_time = [fea_1,fea_2,fea_3,fea_4,fea_5,fea_6];

%% 频域
function fea_spectral = spectral_fea(a)
n = size(a,1);
mag = abs(fft(a));
mag = mag(1:n/2)*2.00/n;

fea_7 = rms(mag); % 均方根频率
fea_8 = std(mag); % 标准差频率
fea_9 = skewness(mag); %偏度
fea_10 = kurtosis(mag); %峰度
fea_11 = mean(power(mag, 3)); %幅值谱功率

fea_spectral = [fea_7,fea_8,fea_9,fea_10,fea_11];

%% 时频
function E_node = energy_fea_com(a) %小波包分量能量
t = wpdec(a,6,'dmey');
%节点数
nodes=get(t,'tn'); %第6层的节点号
N_cfs = length(nodes);
ord=wpfrqord(nodes);  %小波包系数重排，ord是重排后小波包系数索引构成的矩阵　如3层分解的[1;2;4;3;7;8;6;5]
nodes_ord=nodes(ord); %重排后的小波系数
% wavelet packet coefficients. 求取小波包分解的各个节点的小波包系数
for i=1:N_cfs
    cfs6(i,:)=wpcoef(t,nodes_ord(i));
end

for i=1:N_cfs
    E_cfs6(i,:)=norm(cfs6(i,:),2)^2;
end

E_total = sum(E_cfs6);

for i=1:N_cfs
    E_node(i)=100*E_cfs6(i,:)/E_total;
end


function data_fea = extract_fea_ori(data)
num_stat= 11+64;
data_fea = [];
dim_feature = 1;
for i = 1:dim_feature
    data_slice = data;
    data_fea=[time_fea(data_slice),...
        spectral_fea(data_slice),...
        energy_fea_com(data_slice)];
end
data_fea = reshape(data_fea,1, dim_feature*num_stat);


