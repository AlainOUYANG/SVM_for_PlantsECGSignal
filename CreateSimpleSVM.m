%%
% * libsvm toolbox simple usage
%
%% Load data
% * Final data format: m*n, m is number of samples, n is dimension
% * label: m*1 label is -1 and 1
clc
clear
close all
data = load('data_test1.mat');
data = data.data';
% Choose the train number
num_train = 80;
% Create randomly selected sequence
choose = randperm(length(data));
train_data = data(choose(1:num_train),:);
gscatter(train_data(:,1),train_data(:,2),train_data(:,3));
label_train = train_data(:,end);
test_data = data(choose(num_train+1:end),:);
label_test = test_data(:,end);
predict = zeros(length(test_data),1);
%% ----Training model and predict classification
model = svmtrain(label_train,train_data(:,1:end-1),'-t 2');
% -t = 2 Choose radial basis function kernel
true_num = 0;
for i = 1:length(test_data)
    % As the predicton, the first argument of svmpredict can be given randomly
    predict(i) = svmpredict(1,test_data(i,1:end-1),model);
end
%% Show results
figure;
index1 = find(predict==1);
data1 = (test_data(index1,:))';
plot(data1(1,:),data1(2,:),'or');
hold on
index2 = find(predict==-1);
data2 = (test_data(index2,:))';
plot(data2(1,:),data2(2,:),'*');
hold on
indexw = find(predict~=(label_test));
dataw = (test_data(indexw,:))';
plot(dataw(1,:),dataw(2,:),'+g','LineWidth',3);
accuracy = length(find(predict==label_test))/length(test_data);
title(['predict the testing data and the accuracy is :',num2str(accuracy)]);