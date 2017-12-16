%% En-t?te 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DA SILVA BATISTA Rafael Augusto  %
% Signal processing                %
% Last modification : 02/08/17     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
close all;
clear all;

%% Import et traitement du fichier CSV.

% for m = 1 : 1
    
    m=1;
    
    test(m) = m;
    experience = ['chou_notouch_' num2str(test(m)) '.csv'];

    % If you would like to execute for only one file. Uncomment the follow
    % line and comment the main for loop.
    % experience = 't4';

    nom = experience;
    
    % You need to specify a valid adress for get the files.
    nomDuFichierCSV = ['/Users/Ouyangzuokun/Desktop/ProjectEcoGreen/Projet_fichiers_csv/chou_no_touch/' nom];

    fichierCSV = csvread(nomDuFichierCSV, 1, 3);

    fid = fopen(nomDuFichierCSV);
    dateTime = textscan(fid, '%*s %s %*[^\n]','HeaderLines',1, 'Delimiter', ',');
    fclose(fid);

    dateTime = cellstr(dateTime{1,:});
    heureDebut = dateTime{1,1};
    dateTime = datetime(dateTime);
    duration = seconds(dateTime(length(dateTime)) - dateTime(1));

    nombredeMesures = length(fichierCSV(:,1));

    d = fichierCSV(:,7:end);

    dataBrute = [];
    for i = 1:nombredeMesures
        dataBrute = horzcat(dataBrute, d(i,:));
    end
    % save([nom '_signal_brute'], 'dataBrute');
    
    % Frequence d'enchantillonage.
    fs = length(dataBrute)/duration;

    t = 1:(duration-1)/length(dataBrute):duration;
    t = t(1:end-(length(t) - length(dataBrute)));

    %% Analyse spectral du signal et Filtrage
    
    Order = 8;
    
    dataBruteFFT = fft(dataBrute);
    dataBruteFFTConj = dataBruteFFT.*conj(dataBruteFFT)/length(dataBruteFFT);
    dataBruteDB = 20.*log10(dataBruteFFTConj); % Its an amplitude so : AdB = 20.*log10(x)
    f = fs/length(dataBruteFFT)*(0:(length(dataBruteFFT)-1));
    
    %%% Finding la bonne frequence ? analyser.
    
    % From a set value ahead the software will look for the max amplitude in the
    % spectral analyssis. In this case it is set for 3 Hz.
    freqDebutHertz = 3*round(1/(fs/length(dataBruteFFT)));

    % Gettin the max amplitude.
    [maxAmpl, indexMaxAmpl] = max(dataBruteDB(freqDebutHertz:round(length(dataBruteDB)/2)));
    indexMaxAmpl = indexMaxAmpl + freqDebutHertz;
    
    % Vector of 1Hz after the max aplitude
    highLim = dataBruteDB(indexMaxAmpl:indexMaxAmpl+round(1/(fs/length(dataBruteFFT))));
    
    % Getting the high limit of the filter band pass.
    dataBruteDBMeanH = smooth(highLim, round(1/(fs/length(dataBruteFFT))*0.2))'; % 20%
    indexMinAmplLimHigh = length(highLim(dataBruteDBMeanH > -15)); % 15dB = -15dB = 0.17
    if indexMinAmplLimHigh == length(highLim)
        [minAmplLimHigh, indexMinAmplLimHigh] = min(dataBruteDBMeanH);
    end

    % Vector of 1Hz before the max aplitude
    lowLim = dataBruteDB(indexMaxAmpl-round(1/(fs/length(dataBruteFFT))):indexMaxAmpl);
    
    % Getting the low limit of the filter band pass.
    dataBruteDBMeanL = smooth(lowLim, round(1/(fs/length(dataBruteFFT))*0.2))';% 20%
    indexMinAmplLimLow = length(lowLim(dataBruteDBMeanL > -15)); % -15dB = 0.17
    if indexMinAmplLimLow == length(lowLim)
        [minAmplLimLow, indexMinAmplLimLow] = min(dataBruteDBMeanL);
    end

    % Frequency of max amplitude
    fBanPass = f(indexMaxAmpl);
    
    % Getting th cut frequency of the filters.
    fBanPassLimHigh = f(indexMaxAmpl + indexMinAmplLimHigh);
    fBanPassLimLow = f(indexMaxAmpl - (length(dataBruteDBMeanL) - indexMinAmplLimLow));
    
%     % Plotting spectral analysis
%     figure('units','normalized','outerposition',[0 0 1 1])
%     plot(f(indexMaxAmpl - length(lowLim):indexMaxAmpl + length(highLim) - 1), [lowLim highLim])
%     hold on
%     plot(f(indexMaxAmpl - length(lowLim):indexMaxAmpl + length(highLim) - 1), [dataBruteDBMeanL dataBruteDBMeanH])
%     grid on
% 
%     figure('units','normalized','outerposition',[0 0 1 1])
%     subplot(2, 1, 1)
%     plot(f(1:(round(length(f)/10))), dataBruteDB(1:(round(length(dataBruteDB)/10))));
%     title('Analyse spectral de la donn?e brute')
%     xlabel('Fr?quence (Hz)')
%     ylabel('Magnitude (dB)')
%     grid on
%     subplot(2, 1, 2)
%     plot(f(1:(round(length(f)/10))), dataBruteFFTConj(1:(round(length(dataBruteFFTConj)/10))));
%     xlabel('Fr?quence (Hz)')
%     ylabel('Magnitude')
%     grid on
%     % saveas(gcf, [nom '_analyse_spectrale'], 'fig')

    %% Filtering
     
    bpFilt = designfilt('bandpassiir','FilterOrder',Order, ...
        'HalfPowerFrequency1',fBanPassLimLow, ...
        'HalfPowerFrequency2',fBanPassLimHigh, ...
        'SampleRate',fs);

    dataBandFiltre = filter(bpFilt, dataBrute);
    
    lpFilt = designfilt('lowpassiir','FilterOrder',Order, ...
        'PassbandFrequency',fBanPassLimHigh,'PassbandRipple',0.2, ...
        'SampleRate',fs);

     dataLowFiltre = filter(lpFilt, dataBrute);
     
%      save([nom '_signal_filtre_bandPass'], 'dataBandFiltre');
%      save([nom '_signal_filtre_lowPass'], 'dataLowFiltre');
    
    
    %% IDENTIFICATION - BLOCK CODE
    
    % Getting only the postive part of the filtered band-pass signal.
    positifDataFiltre = dataBandFiltre;
    for i = 1 : length(dataBandFiltre)
        if dataBandFiltre(i) < 0
        positifDataFiltre(i) = 0;
        end
    end
    
    % Vector of one second value according to the signal sample frequency.
    unSecond = round(1/((t(end)-t(1))/length(t)));
    
    % Getting the envelope of the filtered band-pass signal using an
    % envelope each 0.5 s.
    
    % In the begin of the signal. There is a issue with the damping of
    % the filter. In this case is set for 4s.
    dampingtTime = 4;
    [envHigh, envLow] = envelope(positifDataFiltre(t > dampingtTime),round(0.5*unSecond),'peak');
    [maxAmplTest, indexMaxAmplTest] = max(envHigh(1:end-dampingtTime*unSecond));
    
    
    % Getting the max of the envelope of the filtered data. 
    % Wich is probably the moment of stress in the plant.
    [envHighTempDebut, envLowTempDebut] = envelope(positifDataFiltre(t <= dampingtTime),round(0.5*unSecond),'peak');
    
    % Temporary vector only for fullfil the begin of the main analysis vector..
    indexAuxTempDebut = length(t(t <= dampingtTime));
    regionOfStressTempDebut = zeros(1, indexAuxTempDebut);
    
    % Getting the stress peak index and relating with time vector.
    indexStressPeak = indexMaxAmplTest;
    tStressPeak = t(indexMaxAmplTest + indexAuxTempDebut);
    
    % Normalising the difference between the max and the mean of the
    % envelope.
    difference = maxAmplTest - mean(envHigh);
    
    %% Mean both sides of the peak for getting the start and the begin of the zone stress.
    
    moyenneApres = sum(envHigh(indexMaxAmplTest:indexMaxAmplTest+unSecond))/length(envHigh(indexMaxAmplTest:indexMaxAmplTest+unSecond));
    moyenneAvant = sum(envHigh(indexMaxAmplTest-unSecond:indexMaxAmplTest))/length(envHigh(indexMaxAmplTest-unSecond:indexMaxAmplTest));
    
    % Getting a parameter normalized for each signal.
    differenceN = moyenneAvant - mean(envHigh);
    
    limLow = unSecond;
    while differenceN > 0.4*difference
        moyenneAvant = sum(envHigh(indexMaxAmplTest-unSecond-limLow:indexMaxAmplTest))/length(envHigh(indexMaxAmplTest-unSecond-limLow:indexMaxAmplTest));
        differenceN = moyenneAvant - mean(envHigh);
        limLow = limLow + unSecond;
    end
    
    % Releasing the difference.
    differenceN = moyenneApres - mean(envHigh);
    
    limHigh = unSecond;
    while differenceN > 0.4*difference
        moyenneApres = sum(envHigh(indexMaxAmplTest:indexMaxAmplTest+unSecond+limHigh))/length(envHigh(indexMaxAmplTest:indexMaxAmplTest+unSecond+limHigh));
        differenceN = moyenneApres - mean(envHigh);
        limHigh = limHigh + unSecond;
    end
    
    % Moment of start and end of the zone stress.
    tDebutStress = t(indexMaxAmplTest + indexAuxTempDebut - limLow);
    tFinStress = t(indexMaxAmplTest + indexAuxTempDebut + limHigh);
    
    % Designing the stress zone vector.
    regionOfStress = ones(1, length(envHigh));
    for i = 1 : length(regionOfStress)
        if indexStressPeak - limLow > i
            regionOfStress(i) = 0;
        end
        if indexStressPeak + limHigh < i
            regionOfStress(i) = 0;
        end
    end

    %% PROCESSING FOR CHARACTERIZE THE STRESS - BLOCK CODE
    
    % Getting the diff of the envelope.
    dataPeak = [ 0 diff(envHigh)];
    
    % Region of stress normalized for the diff.
    regionOfStressDiff = [regionOfStressTempDebut dataPeak].*[regionOfStressTempDebut regionOfStress];
    
    dataPeakTemp = dataPeak.*regionOfStress;
    dataPeakTemp = [regionOfStressTempDebut dataPeakTemp];
    dataPeak = [regionOfStressTempDebut dataPeak];
    
    % Detecting the borders of the stress peak. We get the min and max diff
    % of the stress zone, so it means the "real" start and end of the
    % stress.
    [maxAmplPeak, indexMaxAmplPeak] = max(dataPeakTemp);
    [minAmplPeak, indexMinAmplPeak] = min(dataPeakTemp);
    firstPeakMax = t(indexMaxAmplPeak);
    firstPeakMin = t(indexMinAmplPeak);
    
    % We need the two biggest diffs of the stress zone. So we detect the
    % biggest and then we force it to zero in order to detect the next
    % biggest diff.
    contAuxMaxAmplPeak = indexMaxAmplPeak;
    auxSimetric = 0;
    while dataPeak(contAuxMaxAmplPeak + auxSimetric) >= 0
        dataPeakTemp(contAuxMaxAmplPeak + auxSimetric) = 0;
        dataPeakTemp(contAuxMaxAmplPeak - auxSimetric) = 0;
        auxSimetric = auxSimetric + 1;
    end
    
    contAuxMinAmplPeak = indexMinAmplPeak;
    auxSimetric = 0;
    while dataPeak(contAuxMinAmplPeak + auxSimetric) <= 0
        dataPeakTemp(contAuxMinAmplPeak + auxSimetric) = 0;
        dataPeakTemp(contAuxMinAmplPeak - auxSimetric) = 0;
        auxSimetric = auxSimetric + 1;
    end
    
    [secondMaxAmplPeak, indexSecondMaxAmplPeak] = max(dataPeakTemp);
    [secondMinAmplPeak, indexSecondMinAmplPeak] = min(dataPeakTemp);
    segundoPicoMax = t(indexSecondMaxAmplPeak);
    segundoPicoMin = t(indexSecondMinAmplPeak);
    
    if abs(secondMaxAmplPeak) < 0.5*abs(maxAmplPeak)
        secondMaxAmplPeak = maxAmplPeak;
        indexSecondMaxAmplPeak = indexMaxAmplPeak;
    end
    
    if abs(secondMinAmplPeak) < 0.5*abs(minAmplPeak)
        secondMinAmplPeak = minAmplPeak;
        indexSecondMinAmplPeak = indexMinAmplPeak;
    end
    
    amplitudesVector = [maxAmplPeak minAmplPeak secondMaxAmplPeak secondMinAmplPeak];
    indexAmplitudesVector = [indexMaxAmplPeak indexMinAmplPeak indexSecondMaxAmplPeak indexSecondMinAmplPeak];
    indexAmplitudesVector = indexAmplitudesVector(amplitudesVector ~= 0);
    
    % Vector temporary only for showing a nice view for end user.
    newPeriodStress = zeros(1, length(dataPeak));
    
    maxDataPeak = max(dataPeak);
    maxEnvHigh = max(envHigh);
    minIndexAmplitudesVector = min(indexAmplitudesVector);
    maxIndexAmplitudesVector = max(indexAmplitudesVector);
    
    for j = 1 : length(newPeriodStress)
        if minIndexAmplitudesVector < j && maxIndexAmplitudesVector > j
            newPeriodStress(j) = maxEnvHigh;
        end
    end
    
    % Relating the diff vector with the "real" stress vector. In order to
    % determine where exctly is the stress.
    limNeg = 0;
    limPos = 0;
    for w = 2: length(dataPeak)
        if newPeriodStress(w) == maxDataPeak
            for n = 1 : w-1
                if limPos == 1 && limNeg == 1
                    break
                end
                if w-n > 0
                    if abs(dataPeak(w-n)) > 0.01*maxDataPeak && limNeg == 0
                        newPeriodStress(w) = maxDataPeak;
                        newPeriodStress(w-n) = maxDataPeak;
                    else
                        limNeg = 1;
                    end
                end
                if w+n < length(newPeriodStress)
                    if abs(dataPeak(w+n)) > 0.01*maxDataPeak && limPos == 0
                        newPeriodStress(w) = maxDataPeak;
                        newPeriodStress(w+n) = maxDataPeak;
                    else
                        limPos = 1;
                    end
                end
            end
            limNeg = 0;
            limPos = 0;
        end
    end
    
    % Getting the "real" start and end of the stress. 
    isInverted = 0;
    for p = 2 : length(newPeriodStress)
        if newPeriodStress(p) ~= newPeriodStress(p-1)
            if isInverted == 0
                tStressStart = t(p);
                isInverted = isInverted + 1;
            end
            if isInverted == 1
                tStressEnd = t(p);
            end
        end
    end
    
    tPeriodStress = tStressEnd - tStressStart
   
    if length(newPeriodStress(newPeriodStress ~= 0)) > 5*unSecond
         disp('Flame');
    else
         disp('Cut or touch');
    end
   
    %% Plotting results
     
    regionOfStressNormalise = maxEnvHigh.*regionOfStress;
    
    figure
    plot(t, [envHighTempDebut envHigh], 'b')
    hold on
    plot(t, [regionOfStressTempDebut regionOfStressNormalise], 'r')
    hold on
    plot(t, newPeriodStress, 'g')
    legend('envelope des donn?es', ['debut region stress = ' num2str(tDebutStress) 's' ...
        ', pic du stress = ' num2str(tStressPeak) 'et fin region stress = ' num2str(tFinStress) 's'], ...
        'vrai stress');
    xlabel('Temps (s)')
    ylabel('Amplitude')
    grid on
    
    figure
    plot(t, regionOfStressDiff, 'b')
    legend('Deriv? de l''envelope des donn?es de la r?gion de stress');
    xlabel('Temps (s)')
    ylabel('Amplitude')
    grid on

    figure
    plot(t, dataLowFiltre, 'g');
    hold on
    plot(t, dataBandFiltre, 'b');
    grid on
    title('Signal filtr?')
    xlabel('Temps (s)')
    ylabel('Amplitude (V)')
    legend(['Signal filtr? low-pass ? ' num2str(fBanPassLimHigh) ' Hz'], ...
        ['Signal filtr? band-pass [' num2str(fBanPassLimLow) ' ' ...
        num2str(fBanPassLimHigh) '] Hz' ])

   %% Getting the equation of the stress
    
   zoneEquationStress = zeros(1, length(newPeriodStress));
   
   for j = 1 : length(newPeriodStress)
        if minIndexAmplitudesVector < j && maxIndexAmplitudesVector > j
            zoneEquationStress(j) = 1;
        end
   end
   
   zoneEquationStress = zoneEquationStress.*([envHighTempDebut envHigh]);
   zoneEquationStress = zoneEquationStress(zoneEquationStress ~= 0);
   
   tTest = 1:length(zoneEquationStress);
   
    figure
    plot(tTest, zoneEquationStress, 'r')
    grid on

    FO1 = fit((1:length(zoneEquationStress))', zoneEquationStress', 'poly2');
    FO2 = fit((1:length(zoneEquationStress))', zoneEquationStress', 'poly3');
    FO3 = fit((1:length(zoneEquationStress))', zoneEquationStress', 'poly4');
    FO4 = fit((1:length(zoneEquationStress))', zoneEquationStress', 'poly5');

    hold on
    plot(FO1, 'g--')
    plot(FO2, 'b--')
    plot(FO3, 'm--')
    plot(FO4, 'y--')
    hold off
    legend('data', 'poly2', 'poly3', 'poly4', 'poly5')

%     close all;
% end
% 


