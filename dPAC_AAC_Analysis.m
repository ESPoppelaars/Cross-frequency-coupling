%% Clear workspace and command window
clear
clc

% Start timer
tic

%% Settings

% Sampling rate [Hz]
%InData f = 512;
f = 128;

% Sampling time [s]
T = 1/f;

% Epoch length [s]
EpochLength = 8;

% Time vector [s]
t = 0:T:EpochLength-T;

% The number of samples to cut off from the start and end of every epoch.
% Ideally equal to the filter order of the lower frequency.
edgeLength = 18;

% Number of samples per epoch
epochSamples = EpochLength*f;
epochSamplesCut = epochSamples - 2*edgeLength;

% Number of epochs
epochs = 6;

% Number of samples per file
fileSamples = epochSamples*epochs;
fileSamplesCut = epochSamplesCut*epochs;

% Number of sensors
sensors = 64;

% Number of permutations
permutations = 1000;

% Sensors to combine
combinedSensors = [5 , 38 , 40];
NcombinedSensors = length(combinedSensors);

% Frequency band range [Hz]
%deltaLimits = [1, 4];
%betaLimits = [14, 30];

% Filter settings
DeltaFilter = ButterworthFilter_01_1_4_6_SR_128_24dB;
BetaFilter = ButterworthFilter_12_14_30_32_SR_128_24dB;

%% File processing

% Get csv files to analyse
files = dir('*.csv');

% Initialize results array
% dPAC
export = zeros(length(files),NcombinedSensors);
export_average = zeros(length(files),1);
exportZ = zeros(length(files),NcombinedSensors);
exportZ_average = zeros(length(files),1);
% AAC
exportAACR = zeros(length(files),NcombinedSensors);
exportAACR_average = zeros(length(files),1);
exportAACZ = zeros(length(files),NcombinedSensors);
exportAACZ_average = zeros(length(files),1);
exportAACP = zeros(length(files),NcombinedSensors);
exportAACP_average = zeros(length(files),1);

% Cycle through the files
for fileIndex = 1:length(files)
    
    % Load data into array (skip header)
    data = csvread(files(fileIndex).name,1,0);
    
    % Downsample data from 512 Hz to 128 Hz.
    DownsampledData = downsample(data,4,0);
    
    % Delta band phase variable initialisation
    phaseCombined = zeros(fileSamplesCut*NcombinedSensors,1);
        
    % Delta band amplitude variable initialisation
    amplitudeDeltaCombined = zeros(fileSamplesCut*NcombinedSensors,1);

    % Beta band amplitude variable initialisation
    amplitudeCombined = zeros(fileSamplesCut*NcombinedSensors,1);
    
    % Cycle through sensors
    for sensorIndex = 1:NcombinedSensors
        
        fprintf('Subject %d / %d : Sensor %d / %d\n',fileIndex,length(files),sensorIndex,NcombinedSensors)
        
        % Get the original index of the sensor
        sensor = combinedSensors(sensorIndex);

        % Butterworth IIR bandpass filter. Filter order of 8 for Delta and
        % 34 for Beta. But forwards and backwards to remove phase-shift, 
        % so filter order gets doubled. (16 for delta, 68 for beta.)
        delta = filtfilt(DeltaFilter.SOS,DeltaFilter.ScaleValues,DownsampledData(:,sensor)');
        beta = filtfilt(BetaFilter.SOS,BetaFilter.ScaleValues,DownsampledData(:,sensor)');        
             
        % Delta band phase variable initialisation
        phase = zeros(fileSamplesCut,1);
        amplitudeDelta = zeros(fileSamplesCut,1);

        % Beta band amplitude variable initialisation
        amplitude = zeros(fileSamplesCut,1);
        
        % Cycle through the epochs
        for epoch = 1:epochs
    
            % Start and end indices for the current original epoch
            iStart = (epoch-1)*epochSamples + 1;
            iEnd = epoch*epochSamples;
            
            % Start and end indices for the current cut epoch
            iStartCut = (epoch-1)*epochSamplesCut + 1;
            iEndCut = epoch*epochSamplesCut;
            
            % Delta band phase+ amplitude and beta amplitude through 
            % the Hilbert transform
            phaseEpoch = angle(hilbert(delta(iStart:iEnd)));
            amplitudeDeltaEpoch = abs(hilbert(delta(iStart:iEnd)));
            amplitudeEpoch = abs(hilbert(beta(iStart:iEnd)));

            % Cut the edge effects from the Hilbert-transformed epoch.
            phase(iStartCut:iEndCut) = phaseEpoch(edgeLength + 1:end - edgeLength);
            amplitudeDelta(iStartCut:iEndCut) = amplitudeDeltaEpoch(edgeLength + 1:end - edgeLength);
            amplitude(iStartCut:iEndCut) = amplitudeEpoch(edgeLength + 1:end - edgeLength);
        end
        
        %% AAC
        
        % Get the Pearson correlation coefficient for every sensor and for the 
        % mean of the three sensors.
        R = corr(amplitude,amplitudeDelta);

        % Transform the pearson correlation's r into Fisher's Z (normally
        % distributed).
        Z = 0.5*log((1+R)/(1-R));

        % Calculate two-tailed p-value from z-value, based on a normal
        % distribution.
        P = 2*cdf('Normal',-abs(Z),Z,1);
        
        exportAACR(fileIndex,sensorIndex) = R;
        exportAACZ(fileIndex,sensorIndex) = Z;
        exportAACP(fileIndex,sensorIndex) = P;
        
        %% dPAC
        % Start and end indices for current sensor
        iStart = (sensorIndex-1)*fileSamplesCut + 1;
        iEnd = sensorIndex*fileSamplesCut;
        
        % Append sensors data, to be able to calculate average
        phaseCombined(iStart:iEnd) = phase;
        amplitudeCombined(iStart:iEnd) = amplitude;
        amplitudeDeltaCombined(iStart:iEnd) = amplitudeDelta;
        
        % Phase Clustering bias
        PCbias = mean(exp(1i*phase));

        % Debiased Phase-Ampltiude Cross-Frequency Coupling (dPAC)
        dPAC = mean((exp(1i*phase) - PCbias) .* amplitude);
        
        % Save result in array
        export(fileIndex,sensorIndex) = abs(dPAC);
        
        % Null for z-value
        dPACnull = zeros(1,permutations);

        % Permutate the signal
        for permutation = 1:permutations

            % -- cut-and-paste a random portion of the data; this preserves
            % -- temporal autocorrelation while removing the coupling
            cutLoc = 5 + randperm(length(phase)-10); % -- 5 and 10 prevent the first and last time points from being selected
            cutLoc = cutLoc(1);
            phaseShuffled = phase([cutLoc:end 1:cutLoc-1]);

            % Compute surrogate dPAC
            dPACnull(permutation) = abs(mean((exp(1i*phaseShuffled) - mean(exp(1i*phaseShuffled))) .* amplitude));
        end

        % dPAC z-score
        meandPACnull = mean(dPACnull);
        stddPACnull = std(dPACnull);
        dPACz = (abs(dPAC) - meandPACnull) ./ stddPACnull;
        
        % Save result in array
        exportZ(fileIndex,sensorIndex) = abs(dPACz);
    end
end

    %% dPAC for the average of the sensors
    
    % Average the dPAC and Z-value per participant
    export_average = mean(export,2);
    exportZ_average = mean(exportZ,2);
    
    %% AAC for the average of the sensors
    
    % Average the Pearson correlation coefficient for every sensor combined.
    exportAACR_average = mean(exportAACR,2);

    for fileindex = 1:length(files)
        
        R = exportAACR_average(fileindex);
        
        % Transform the pearson correlation's r into Fisher's Z (normally
        % distributed).
        Z = 0.5*log((1+R)/(1-R));

        % Calculate two-tailed p-value from z-value, based on a normal
        % distribution.
        P = 2*cdf('Normal',-abs(Z),Z,1);

        exportAACZ_average(fileIndex) = Z;
        exportAACP_average(fileIndex) = P;
    end

%% Export results

% Write CSV
csvwrite('000_IndividualResults_dPAC_128Hz.txt',[export,export_average]);
csvwrite('000_IndividualResults_z_128Hz.txt',[exportZ,exportZ_average]);
csvwrite('000_IndividualResults_AAC_R_128Hz.txt',[exportAACR,exportAACR_average]);
csvwrite('000_IndividualResults_AAC_P_128Hz.txt',[exportAACP,exportAACP_average]);
csvwrite('000_IndividualResults_AAC_Z_128Hz.txt',[exportAACZ,exportAACZ_average]);

%% Polar plot of phase angles

% Downsample phase into 1/400th its original size, to plot.
downSampleRatio = round(length(phaseCombined)/300);

phasePlot = downsample(phaseCombined,downSampleRatio);
amplitudePlot = downsample(amplitudeCombined,downSampleRatio);

Colour = [.6,.6,.6];
LineWidth = .7;
ampZoom = 10;

figure(1)
clf
polarplot([0,phasePlot(1)],[0,min(amplitudePlot(1),ampZoom)] , 'Color',Colour , 'LineWidth',LineWidth)
hold on
for i = 2:length(phasePlot)
    polarplot([0,phasePlot(i)],[0,min(amplitudePlot(i),ampZoom)] , 'Color',Colour , 'LineWidth',LineWidth)
end
polarplot([0,mean(angle(dPAC))],[0,mean(abs(dPAC))],'k' , 'LineWidth',4)

saveas(figure(1),'PolarPlot','png')
hold off

%%
% Stop timer
toc

% Notify user with sound
beep