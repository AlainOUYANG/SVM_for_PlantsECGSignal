figure;
load('/Users/Ouyangzuokun/Documents/MATLAB/Projet_EchoGreen/RawData/chou_notouch/chou_notouch_11_raw_signal.mat');
subplot(2,2,1); plot(0:18/length(rawData):18-18/length(rawData),rawData);
title('Raw data for chou\_notouch\_11'); xlabel('Time (s)'); ylabel('Amplitude (V)');

load('/Users/Ouyangzuokun/Documents/MATLAB/Projet_EchoGreen/RawData/chou_longtouch/chou_longtouch_38_raw_signal.mat');
subplot(2,2,2); plot(0:18/length(rawData):18-18/length(rawData),rawData);
title('Raw data for chou\_longtouch\_38'); xlabel('Time (s)'); ylabel('Amplitude (V)');

load('/Users/Ouyangzuokun/Documents/MATLAB/Projet_EchoGreen/RawData/chou_shorttouch/chou_shorttouch_25_raw_signal.mat');
subplot(2,2,3); plot(0:18/length(rawData):18-18/length(rawData),rawData);
title('Raw data for chou\_shorttouch\_25'); xlabel('Time (s)'); ylabel('Amplitude (V)');

load('/Users/Ouyangzuokun/Documents/MATLAB/Projet_EchoGreen/RawData/chou_water/chou_water_66_raw_signal.mat');
subplot(2,2,4); plot(0:18/length(rawData):18-18/length(rawData),rawData);
title('Raw data for chou\_water\_66'); xlabel('Time (s)'); ylabel('Amplitude (V)');