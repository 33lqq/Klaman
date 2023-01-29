data=[A;B;C;D;E;F;G;H];
data_nominal = data;
% data_nominal = mapminmax(data);
%% Label
label(1:15000,:) = ones(15000,1)*1;
label(15000+1:15000*2,:) = ones(15000,1)*2;
label(15000*2+1:15000*3,:)=ones(15000,1)*3;
label(15000*3+1:15000*4,:)=ones(15000,1)*4;
label(15000*4+1:15000*5,:)=ones(15000,1)*5;
label(15000*5+1:15000*6,:)=ones(15000,1)*6;
label(15000*6+1:15000*7,:)=ones(15000,1)*7;
label(15000*7+1:15000*8,:)=ones(15000,1)*8;

data_nominal = [data_nominal(:,1),data_nominal(:,2),data_nominal(:,3),data_nominal(:,6),data_nominal(:,9),data_nominal(:,11+1),data_nominal(:,11+2),data_nominal(:,11+3),data_nominal(:,11+21),data_nominal(:,11+26),data_nominal(:,11+33),data_nominal(:,11+44),data_nominal(:,11+49),data_nominal(:,11+64),...
data_nominal(:,(11+64)*1+2),data_nominal(:,(11+64)*1+3),data_nominal(:,(11+64)*1+11+14),data_nominal(:,(11+64)*1+11+17),data_nominal(:,(11+64)*1+11+30),data_nominal(:,(11+64)*1+11+31),data_nominal(:,(11+64)*1+11+38),data_nominal(:,(11+64)*1+11+45),data_nominal(:,(11+64)*1+11+51),data_nominal(:,(11+64)*1+11+60),...
data_nominal(:,(11+64)*2+1),data_nominal(:,(11+64)*2+2),data_nominal(:,(11+64)*2+3),data_nominal(:,(11+64)*2+6),data_nominal(:,(11+64)*2+9),data_nominal(:,(11+64)*2+11+8),data_nominal(:,(11+64)*2+11+17),data_nominal(:,(11+64)*2+11+18),data_nominal(:,(11+64)*2+11+21),data_nominal(:,(11+64)*2+11+25),data_nominal(:,(11+64)*2+11+38),data_nominal(:,(11+64)*2+11+39),data_nominal(:,(11+64)*2+11+53),data_nominal(:,(11+64)*2+11+56),data_nominal(:,(11+64)*2+11+58)];
data_total = [data_nominal,label];
rowrank = randperm(size(data_total, 1)); % size获得a的行数，randperm打乱各行的顺序
data_total = data_total(rowrank,:);              % 按照rowrank重新排列各行，注意rowrank的位置

