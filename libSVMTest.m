clc; clear; close all;

data = [1.5 1234; 1.3 1024; 1.9 8765; 2.2 8989];
label = [1;1;-1;-1];

model = svmtrain(label, data);

model
Parameters = model.Parameters
Label = model.Label
nr_class = model.nr_class
totalSV = model.totalSV
nSV = model.nSV 

testdata = [2.7034e-05 8922];
testdatalabel = 1;

[predictlabel] = svmpredict(testdatalabel, testdata, model);

predictlabel