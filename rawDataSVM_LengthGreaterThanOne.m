tic; clc; clear; close all;

%% Calculate the max of each rawData and save as the dataset with corresopding labels

dataNum = 200;
data = zeros(dataNum, 2);
label = ones(dataNum, 1); % label = 1, it's long touch; label = -1, it's short touch

for i=1:dataNum
    if i <= 100
        % Load rawData
        test(i) = i;
        dataName = ['/Users/Ouyangzuokun/Documents/MATLAB/Projet_EchoGreen/RawData/chou_longtouch/chou_longtouch_' num2str(test(i)) '_raw_signal.mat'];
        load(dataName);

        % %calculate the maximum of each rawData, if max > 1, label = 1, else, label = -1
        data(i, 1) = max(rawData);
        data(i, 2) = length(find(rawData>1));
        label(i) = 1;
    else
        test(i-100) = i-100;
        dataName = ['/Users/Ouyangzuokun/Documents/MATLAB/Projet_EchoGreen/RawData/chou_shorttouch/chou_shorttouch_' num2str(test(i-100)) '_raw_signal.mat'];
        load(dataName);

        % %calculate the maximum of each rawData, if max > 1, label = 1, else, label = -1
        data(i, 1) = max(rawData);
        data(i, 2) = length(find(rawData>1));
        label(i) = -1;
    end
end

%% svm train
model = svmtrain(label, data, '-s 0 -t 2 -c 1.2 -g 2.8');

model
% Parameters = model.Parameters
% Label = model.Label
% nr_class = model.nr_class
% totalSV = model.totalSV
% nSV = model.nSV 

%% Use trained model to see its effect on the classification of the training set
[ptrain] = svmpredict(label, data, model);

%% Predict test data
% testdata = 1.231364;
% load('/Users/Ouyangzuokun/Documents/MATLAB/Projet_EchoGreen/RawData/chou_longtouch/chou_longtouch_93_raw_signal.mat');
% 
% testdata = [max(rawData) length(find(rawData>1))];
% testlabel = 1;
% 
% [ptest] = svmpredict(testlabel, testdata, model);
% 
% ptest


%% Load multiple test data
testNum = 100;
testdata = zeros(testNum, 2);
testlabel = ones(testNum, 1); % label = 1, there is some stress; label = -1, there is no stress

for i=1:testNum
    if i <= 50
        % Load rawData
        test(i) = i;
        dataName = ['/Users/Ouyangzuokun/Documents/MATLAB/Projet_EchoGreen/RawData/chou_longtouch/chou_longtouch_' num2str(test(i)) '_raw_signal.mat'];
        load(dataName);

        % calculate the maximum of each rawData, if max > 1, label = 1, else, label = -1
        testdata(i, 1) = max(rawData);
        testdata(i, 2) = length(find(rawData>1));

        testlabel(i) = 1;
    else
        test(i+50) = i+50;
        dataName = ['/Users/Ouyangzuokun/Documents/MATLAB/Projet_EchoGreen/RawData/chou_shorttouch/chou_shorttouch_' num2str(test(i+50)) '_raw_signal.mat'];
        load(dataName);

        % calculate the maximum of each rawData, if max > 1, label = 1, else, label = -1
        testdata(i, 1) = max(rawData);
        testdata(i, 2) = length(find(rawData>1));

        testlabel(i) = -1;
    end
end

[ptest] = svmpredict(testlabel, testdata, model);

% ptest
%%
toc;



