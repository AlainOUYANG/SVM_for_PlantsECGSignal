%% Get raw data from csv files

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Zuokun OUYANG         %
% Last midification: 04/12/2017 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear;
close all;

%% Import csv file and export raw signal data in matrix
FileNum = 150;
FileNameStart = 'chou_shorttouch_';
FilePath = '/Users/Ouyangzuokun/Desktop/ProjectEcoGreen/Projet_fichiers_csv/chou_shorttouch/';
for m = 1:FileNum
    
    test(m) = m;
    csvFileName = [FileNameStart num2str(test(m)) '.csv'];
    
    FullCSVName = [FilePath csvFileName];
    
    fileCSV = csvread(FullCSVName, 1, 3);
    
    fid = fopen(FullCSVName);
    dateTime = textscan(fid, '%*s %s %*[^\n]','HeaderLines',1, 'Delimiter', ',');
    fclose(fid);
    
    dateTime = cellstr(dateTime{1,:});
    startTime = dateTime{1,1};
    dateTime = datetime(dateTime);
    duration = seconds(dateTime(length(dateTime)) - dateTime(1));
    
    MesuresNum = length(fileCSV(:,1));
    
    d = fileCSV(:,7:end);
    
    rawData = [];
    for i = 1:MesuresNum
        rawData = horzcat(rawData, d(i,:));
    end
    save([FileNameStart num2str(test(m)) '_raw_signal'], 'rawData');
end