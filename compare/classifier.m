N =200;
label(1:N,:) = ones(N,1)*1;
label(N+1:2*N,:) = ones(N,1)*2;
label(2*N+1:3*N,:)=ones(N,1)*3;
label(3*N+1:4*N,:)=ones(N,1)*4;
label(4*N+1:5*N,:)=ones(N,1)*5;
label(5*N+1:6*N,:)=ones(N,1)*6;
label(6*N+1:7*N,:)=ones(N,1)*7;
label(7*N+1:8*N,:)=ones(N,1)*8;

data=[A;B;C;D;E;F;G;H];

data_nominal = mapminmax(data);
data_nominal = [data_nominal(:,1),data_nominal(:,2),data_nominal(:,3),data_nominal(:,6),data_nominal(:,9),data_nominal(:,11+1),data_nominal(:,11+2),data_nominal(:,11+3),data_nominal(:,11+21),data_nominal(:,11+26),data_nominal(:,11+33),data_nominal(:,11+44),data_nominal(:,11+49),data_nominal(:,11+64),...
    data_nominal(:,(11+64)*1+2),data_nominal(:,(11+64)*1+3),data_nominal(:,(11+64)*1+11+14),data_nominal(:,(11+64)*1+11+17),data_nominal(:,(11+64)*1+11+30),data_nominal(:,(11+64)*1+11+31),data_nominal(:,(11+64)*1+11+38),data_nominal(:,(11+64)*1+11+45),data_nominal(:,(11+64)*1+11+51),data_nominal(:,(11+64)*1+11+60),...
    data_nominal(:,(11+64)*2+1),data_nominal(:,(11+64)*2+2),data_nominal(:,(11+64)*2+3),data_nominal(:,(11+64)*2+6),data_nominal(:,(11+64)*2+9),data_nominal(:,(11+64)*2+11+8),data_nominal(:,(11+64)*2+11+17),data_nominal(:,(11+64)*2+11+18),data_nominal(:,(11+64)*2+11+21),data_nominal(:,(11+64)*2+11+25),data_nominal(:,(11+64)*2+11+38),data_nominal(:,(11+64)*2+11+39),data_nominal(:,(11+64)*2+11+53),data_nominal(:,(11+64)*2+11+56),data_nominal(:,(11+64)*2+11+58)];

data_total = [data_nominal,label];

rowrank = randperm(size(data_total, 1)); % size获得a的行数，randperm打乱各行的顺序
data_total = data_total(rowrank,:);              % 按照rowrank重新排列各行，注意rowrank的位置