clc; clear; close all;
%%
load('/Users/Ouyangzuokun/Documents/MATLAB/Projet_EchoGreen/RawData/chou_notouch/chou_notouch_12_raw_signal.mat');
figure;plot(abs(fftshift(fft(rawData))));
temp = find(rawData>1);
fprintf('notouch: length of stress: %d\n', length(temp))

%%
load('/Users/Ouyangzuokun/Documents/MATLAB/Projet_EchoGreen/RawData/chou_shorttouch/chou_shorttouch_80_raw_signal.mat');
figure;plot(abs(fftshift(fft(rawData))));
temp = find(rawData>1);
fprintf('shorttouch: length of stress: %d, max interval: %d\n', length(temp), temp(length(temp)) - temp(1))

%%
load('/Users/Ouyangzuokun/Documents/MATLAB/Projet_EchoGreen/RawData/chou_longtouch/chou_longtouch_17_raw_signal.mat');
figure;plot(abs(fftshift(fft(rawData))));
temp = find(rawData>1);
fprintf('longtouch: length of stress: %d, max interval: %d\n', length(temp), temp(length(temp)) - temp(1))

%%
load('/Users/Ouyangzuokun/Documents/MATLAB/Projet_EchoGreen/RawData/chou_water/chou_water_39_raw_signal.mat');
figure;plot(abs(fftshift(fft(rawData))));
temp = find(rawData>1);
fprintf('water: length of stress: %d, max interval: %d\n', length(temp), temp(length(temp)) - temp(1))