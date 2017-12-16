%% 
% * Logistic method to Regression Analysis and Classification Design
% * Multi-classes classification
% 
%% Create 200 non-linear datas, 2 types
clc
clear
close all
data1 = rand(1,1000)*4 + 1;
data2 = rand(1,1000)*4 + 1;
circle_inx = data1 - 3;
circle_iny = data2 - 3;
r_in = circle_inx.^2 + circle_iny.^2;
index_in = find(r_in<0.8);
data_in = [data1(index_in(1:100));data2(index_in(1:100));-1*ones(1,100)];

data1 = rand(1,1000)*4 + 1;
data2 = rand(1,1000)*4 + 1;
circle_inx = data1 - 3;
circle_iny = data2 - 3;
r_in = circle_inx.^2 + circle_iny.^2;
index_out = find((r_in>1.2)&(r_in<4));
data_out = [data1(index_out(1:100));data2(index_out(1:100));ones(1,100)];
plot(data_in(1,:),data_in(2,:),'or');
hold on;
plot(data_out(1,:),data_out(2,:),'*');
data = [data_in,data_out];
save data_test1.mat data