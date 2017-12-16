%% libsvm tutorial
tic;
close all;
clear;
clc;
format compact;

%%
% Load data
[heart_scale_label, heart_scale_inst] = libsvmread('heart_scale');
data = heart_scale_inst;
label = heart_scale_label;

%%
% Use the first 200 data for training and the last 70 data for testing
ind = 200;
traindata = data(1:ind, :);
trainlabel = label(1:ind, :);
testdata = data(ind+1:end, :);
testlabel = label(ind+1:end, :);

% Use the trian data to create the classifier model
model = svmtrain(trainlabel, traindata, '-s 0 -t 2 -c 1.2 -g 2.8');

% Classifier model
model
Parameters = model.Parameters
Label = model.Label
nr_class = model.nr_class
totalSV = model.totalSV
nSV = model.nSV

%% Use trained model to see its effect on the classification of the training set
[ptrain] = svmpredict(trainlabel, traindata, model);

%% Predict test data
[ptest] = svmpredict(testlabel, testdata, model);

%%
toc;