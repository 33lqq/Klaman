% y=C*xupd;
% y=C*xupd_1;
% y=C*xupd_2;
% xupd_3 = reshape(xupd_3,4,1501);
% y=C*xupd_3;
% xupd_4 = reshape(xupd_4,4,1501);
% y=C*xupd_4;
y=Ytrue;
for i=1:30
    %% 1
    data1 = y(1,2+(i-1)*50:i*50+1);
    X(1) = rms(data1);
    X(2) = var(data1);
    X(3) = max(data1);
     X(4) = max(data1) - min(data1);

     n = size(data1,2);
%     mag = abs(fft(data1)./n*2);
    mag = abs(fft(data1));

    mag = mag(1:n/2)*2.00/n;
    X(5) = skewness(mag); %偏度
    t = wpdec(data1,6,'dmey');
    %节点数
    nodes=get(t,'tn'); %第6层的节点号
    N_cfs = length(nodes);
    ord=wpfrqord(nodes);  %小波包系数重排，ord是重排后小波包系数索引构成的矩阵　如3层分解的[1;2;4;3;7;8;6;5]
    nodes_ord=nodes(ord); %重排后的小波系数
    % wavelet packet coefficients. 求取小波包分解的各个节点的小波包系数
    for j=1:N_cfs
        cfs6(j,:)=wpcoef(t,nodes_ord(j));
    end
    
    for j=1:N_cfs
        E_cfs6(j,:)=norm(cfs6(j,:),2)^2;
    end
    
    E_total = sum(E_cfs6);
    
    for j=1:N_cfs
        E_node(j)=100*E_cfs6(j,:)/E_total;
    end
    X(6) = E_node(1);
    X(7) = E_node(2);
    X(8)= E_node(3);
    X(9)= E_node(21);
    X(10)= E_node(26);
    X(11)= E_node(33);
    X(12)= E_node(44);
    X(13)= E_node(49);
    X(14)= E_node(64);
%% 2
    data2 = y(2,2+(i-1)*50:i*50+1);

    X(15) = var(data2);
    X(16) = max(data2);
    t = wpdec(data2,6,'dmey');
    %节点数
    nodes=get(t,'tn'); %第6层的节点号
    N_cfs = length(nodes);
    ord=wpfrqord(nodes);  %小波包系数重排，ord是重排后小波包系数索引构成的矩阵　如3层分解的[1;2;4;3;7;8;6;5]
    nodes_ord=nodes(ord); %重排后的小波系数
    % wavelet packet coefficients. 求取小波包分解的各个节点的小波包系数
    for j=1:N_cfs
        cfs6(j,:)=wpcoef(t,nodes_ord(j));
    end
    
    for j=1:N_cfs
        E_cfs6(j,:)=norm(cfs6(j,:),2)^2;
    end
    
    E_total = sum(E_cfs6);
    
    for j=1:N_cfs
        E_node(j)=100*E_cfs6(j,:)/E_total;
    end
    X(17) = E_node(14);
    X(18) = E_node(17);
    X(19)= E_node(30);
    X(20)= E_node(31);
    X(21)= E_node(38);
    X(22)= E_node(45);
    X(23)= E_node(51);
    X(24)= E_node(60);
%% 3
    data3 = y(3,2+(i-1)*50:i*50+1);
    X(25) = rms(data3);
    X(26) = var(data3);
    X(27) = max(data3);
    X(28) = max(data3) - min(data3);
    n = size(data3,2);
%     mag = abs(fft(data3)./n*2);
     mag = abs(fft(data3));

    mag = mag(1:n/2)*2.00/n;
    X(29) = skewness(mag); %偏度
    t = wpdec(data3,6,'dmey');
    %节点数
    nodes=get(t,'tn'); %第6层的节点号
    N_cfs = length(nodes);
    ord=wpfrqord(nodes);  %小波包系数重排，ord是重排后小波包系数索引构成的矩阵　如3层分解的[1;2;4;3;7;8;6;5]
    nodes_ord=nodes(ord); %重排后的小波系数
    % wavelet packet coefficients. 求取小波包分解的各个节点的小波包系数
    for j=1:N_cfs
        cfs6(j,:)=wpcoef(t,nodes_ord(j));
    end
    
    for j=1:N_cfs
        E_cfs6(j,:)=norm(cfs6(j,:),2)^2;
    end
    
    E_total = sum(E_cfs6);
    
    for j=1:N_cfs
        E_node(j)=100*E_cfs6(j,:)/E_total;
    end
    X(30) = E_node(8);
    X(31) = E_node(17);
    X(32) = E_node(18);
    X(33)= E_node(21);
    X(34)= E_node(25);
    X(35)= E_node(38);
    X(36)= E_node(39);
    X(37)= E_node(53);
    X(38)= E_node(56);
    X(39)= E_node(58);
    
%     X = mapminmax(X);
%     yfit(i) = MiddleNN.predictFcn(X);
    yfit(i) = QuadSVM.predictFcn(X);
%     yfit(i) = trainedModel.predictFcn(X);
end