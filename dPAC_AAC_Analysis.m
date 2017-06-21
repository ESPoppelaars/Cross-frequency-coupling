%% Recommendations:
% Needs the toolboxes: parallel computing and signal processing.

%% Clear workspace and command window
clear
clc

% Start timer
tic

%% Settings

% Sampling rate [Hz]
%In data: 512 Hz, but downsampled in code here;
f = 128;

% Sampling time [s]
T = 1/f;

% Time vector [s]
t = 0:T:8-T;

% The number of samples to cut off from the start and end of every epoch.
% Ideally equal to the filter order of the lower frequency.
edgeLength = 16;

% Number of samples per epoch (seconds*frequency)
epochSamples = 8*f;
epochSamplesCut = epochSamples - 2*edgeLength;

% Number of epochs.
epochs = 6;

% Number of samples per file
fileSamples = epochSamples*epochs;
fileSamplesCut = epochSamplesCut*epochs;

% Number of sensors
sensors = 64;

% Number of permutations dPAC
permutations = 1000;

% Sensors to combine
combinedSensors = [5 , 38 , 40];
NcombinedSensors = length(combinedSensors);

% Filters designed and exported with fdatool: 
% delta = 1 - 4 Hz, beta = 14 - 30 Hz
DeltaFilter = DeltaFilterDesign;
BetaFilter = BetaFilterDesign;

%% File processing

% Get csv files from current directory
files = dir('*.csv');

% Number of files / subjects
Nfiles = length(files);

% Initialise matrices for appending all phase and amplitude measures
phaseAppended = zeros(Nfiles*fileSamplesCut,NcombinedSensors);
amplitudeAppended = zeros(Nfiles*fileSamplesCut,NcombinedSensors);
amplitudeDeltaAppended = zeros(Nfiles*fileSamplesCut,NcombinedSensors);
AAC_corr_appended = zeros(Nfiles*epochs,NcombinedSensors);

% Cycle through the files
for fileIndex = 1:Nfiles
    
    % Show progress message
    fprintf('Loading subject %d / %d\n',fileIndex,length(files))
    
    % Load data into array (skip header)
    data = csvread(files(fileIndex).name,1,0);
    
    % Downsample data from 512 Hz to 128 Hz.
    DownsampledData = downsample(data,4,0);
    
    % Cycle through sensors
    for sensorIndex = 1:NcombinedSensors
        
        % Get the original index of the sensor
        sensor = combinedSensors(sensorIndex);
  
        % Butterworth IIR bandpass filter. Filter order of 8 for Delta and
        % 34 for Beta. But forwards and backwards to remove phase-shift, 
        % so filter order gets doubled: 16 for delta, 68 for beta.
        delta = filtfilt(DeltaFilter.SOS,DeltaFilter.ScaleValues,DownsampledData(:,sensor)');
        beta = filtfilt(BetaFilter.SOS,BetaFilter.ScaleValues,DownsampledData(:,sensor)');        

        % Delta band phase variable initialisation
        phase  = zeros(fileSamplesCut,1);
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
            
            amplitudeDelta(iStartCut:iEndCut) = (amplitudeDeltaEpoch(edgeLength + 1:end - edgeLength))/max(amplitudeDeltaEpoch(edgeLength + 1:end - edgeLength));
            amplitude(iStartCut:iEndCut) = (amplitudeEpoch(edgeLength + 1:end - edgeLength))/max(amplitudeEpoch(edgeLength + 1:end - edgeLength));
                        
            % Calculate the correlation between the mean envelope of beta
            % and the amplitude of delta.
            AAC_corr_current = corr(amplitude(iStartCut:iEndCut),amplitudeDelta(iStartCut:iEndCut));
            AAC_corr(epoch,1) = AAC_corr_current;
                      
        end
              
        % Start and end indices for current sensor (used for appending).
        iStartCut = (fileIndex-1)*fileSamplesCut + 1;
        iEndCut = fileIndex*fileSamplesCut;
        iStartCut_AAC = (fileIndex-1)*epochs + 1;
        iEndCut_AAC = fileIndex*epochs;
        
        % Append sensors data, to be able to calculate average
        phaseAppended(iStartCut:iEndCut,sensorIndex) = phase;
        amplitudeDeltaAppended(iStartCut:iEndCut,sensorIndex) = amplitudeDelta;
        amplitudeAppended(iStartCut:iEndCut,sensorIndex) = amplitude; 
        AAC_corr_appended(iStartCut_AAC:iEndCut_AAC,sensorIndex) = AAC_corr;
                
    end
end

%Create a combined variable of all three sensors for the amplitudes of
%Delta and Beta.
amplitudeDeltaAppendedAll = [amplitudeDeltaAppended(:,1) ; amplitudeDeltaAppended(:,2); amplitudeDeltaAppended(:,3)];
amplitudeAppendedAll = [amplitudeAppended(:,1) ; amplitudeAppended(:,2); amplitudeAppended(:,3)];

% Append all correlations of each sensor into a combined variable.
AAC_corr_Comb = [AAC_corr_appended(:,1) ; AAC_corr_appended(:,2) ; AAC_corr_appended(:,3)];

% Calculate the standard deviation and standard error of the mean.
AAC_corr_SD_Comb = [std(AAC_corr_appended(:,1)) ; std(AAC_corr_appended(:,2)) ; std(AAC_corr_appended(:,3))];
AAC_SD = mean(AAC_corr_SD_Comb);
AAC_SEM = AAC_SD / sqrt(Nfiles);

%% Calculate significance of group amplitude-amplitude coupling (according to Knyazev 2011)

% Transform the pearson correlation's r into Fisher's Z (normally
% % distributed) for each epoch and participant.
Z1_separate = 0.5.*log((1+AAC_corr_appended(:,1))./(1-AAC_corr_appended(:,1)));
Z2_separate = 0.5.*log((1+AAC_corr_appended(:,2))./(1-AAC_corr_appended(:,2)));
Z3_separate = 0.5.*log((1+AAC_corr_appended(:,3))./(1-AAC_corr_appended(:,3)));
ZM_separate = 0.5.*log((1+AAC_corr_Comb)./(1-AAC_corr_Comb));

% Calculate the mean Fisher Z.
Z1_M_separate = mean(Z1_separate);
Z2_M_separate = mean(Z2_separate);
Z3_M_separate = mean(Z3_separate);
ZM_M_separate = mean(ZM_separate);

% Tranform back to the mean correlation coefficient based on the mean
% Fisher Z.
R1_separate = (exp(2*Z1_M_separate)-1)/(exp(2*Z1_M_separate)+1);
R2_separate = (exp(2*Z2_M_separate)-1)/(exp(2*Z2_M_separate)+1);
R3_separate = (exp(2*Z3_M_separate)-1)/(exp(2*Z3_M_separate)+1);
RM_separate = (exp(2*ZM_M_separate)-1)/(exp(2*ZM_M_separate)+1);

% Calculate two-tailed p-value from z-value, based on a normal
% distribution.
P1_separate = 2*cdf('Normal',-abs(Z1_M_separate),Z1_M_separate,1/(sqrt(Nfiles-3)));
P2_separate = 2*cdf('Normal',-abs(Z2_M_separate),Z2_M_separate,1/(sqrt(Nfiles-3)));
P3_separate = 2*cdf('Normal',-abs(Z3_M_separate),Z3_M_separate,1/(sqrt(Nfiles-3)));
PM_separate = 2*cdf('Normal',-abs(ZM_M_separate),ZM_M_separate,1/(sqrt(Nfiles-3)));

%% Individual AAC.

% Tranform back each correlation coefficient based on each Fisher Z.
R_PerSubjectAndEpoch1 = (exp(2.*Z1_separate)-1)./(exp(2.*Z1_separate)+1);
R_PerSubjectAndEpoch2 = (exp(2.*Z2_separate)-1)./(exp(2.*Z2_separate)+1);
R_PerSubjectAndEpoch3 = (exp(2.*Z3_separate)-1)./(exp(2.*Z3_separate)+1);

% Calculate the mean correlation coeffcient per participant for each sensor separately.
for subject = 1:Nfiles
    R_PerSubject1(subject) = mean(R_PerSubjectAndEpoch1(subject*1:subject*epochs,1));
    R_PerSubject2(subject) = mean(R_PerSubjectAndEpoch2(subject*1:subject*epochs,1));
    R_PerSubject3(subject) = mean(R_PerSubjectAndEpoch3(subject*1:subject*epochs,1));
end

% Append the (transposed) sensors and calculate the mean correlation coefficient for each participant.
R_PerSubjectAll = [R_PerSubject1',R_PerSubject2',R_PerSubject3'];
R_PerSubjectMean = mean(R_PerSubjectAll,2);

%% Export AAC results.

%Write CSV textfiles with AAC results (4 numbers, for sensors 1-3 and the
%combined one), one document for r coefficients, p values, and z-transformed values.
%Group.
csvwrite('000_GroupResults_AAC_R.txt',[R1_separate,R2_separate,R3_separate,RM_separate]);
csvwrite('000_GroupResults_AAC_P.txt',[P1_separate,P2_separate,P3_separate,PM_separate]);
csvwrite('000_GroupResults_AAC_Z.txt',[Z1_M_separate,Z2_M_separate,Z3_M_separate,ZM_M_separate]);
csvwrite('000_GroupResults_AAC_SD_SEM.txt',[AAC_SD,AAC_SEM]);
%Individual.
csvwrite('000_GroupResults_AAC_R_PerSubject.txt',[R_PerSubjectMean]);

%% Calculate dPAC results (according to Van Driel 2015)

% Initialise results matrices
dPACresults = zeros(1,4);
Zresults = zeros(1,4);

% Initialise matrices for appending all phase and amplitude measures
phaseCombined = zeros(Nfiles*fileSamplesCut*NcombinedSensors,1);
amplitudeCombined = zeros(Nfiles*fileSamplesCut*NcombinedSensors,1);

% Cycle through the sensors
for sensor = 1:NcombinedSensors
    
    % Print empty line
    fprintf('\n');
    
    % Start and end indices for current sensor
    iStart = (sensor-1)*Nfiles*fileSamplesCut + 1;
    iEnd = sensor*Nfiles*fileSamplesCut;

    % Append sensors data, to be able to calculate average
    phaseCombined(iStart:iEnd) = phaseAppended(:,sensor);
    amplitudeCombined(iStart:iEnd) = amplitudeAppended(:,sensor); 
    
    % Phase Clustering bias
    PCbias = mean(exp(1i*phaseAppended(:,sensor)));

    % Debiased Phase-Ampltiude Cross-Frequency Coupling (dPAC)
    dPAC = mean((exp(1i*phaseAppended(:,sensor)) - PCbias) .* amplitudeAppended(:,sensor));
    
    % Store result
    dPACresults(sensor) = abs(dPAC);
    
    % Null for z-value
    dPACnull = zeros(1,permutations);

    % Permutate the signal
    for permutation = 1:permutations
        
        % Show progress
        if mod(permutation,100) == 0
            fprintf('Sensor %d / %d : Permutation %d / %d\n',sensor,NcombinedSensors,permutation,permutations)
        end

        % -- cut-and-paste a random portion of the data; this preserves
        % -- temporal autocorrelation while removing the coupling
        cutLoc = 5 + randperm(length(phaseAppended(:,sensor))-10); % -- 5 and 10 prevent the first and last time points from being selected
        cutLoc = cutLoc(1);
        phaseShuffled = phaseAppended([cutLoc:end 1:cutLoc-1],sensor);

        % Compute surrogate dPAC
        dPACnull(permutation) = abs(mean((exp(1i*phaseShuffled) - mean(exp(1i*phaseShuffled))) .* amplitudeAppended(:,sensor)));
    end

    % dPAC z-score
    meandPACnull = mean(dPACnull);
    stddPACnull = std(dPACnull);
    dPACz = (abs(dPAC) - meandPACnull) ./ stddPACnull;
    
    % Store result
    Zresults(sensor) = dPACz;
end


%% Combined results

% Print empty line
fprintf('\n');

% Phase Clustering bias
PCbias = mean(exp(1i*phaseCombined));

% Debiased Phase-Amplitude Cross-Frequency Coupling (dPAC)
dPAC = mean((exp(1i*phaseCombined) - PCbias) .* amplitudeCombined);

% Store result
dPACresults(NcombinedSensors + 1) = abs(dPAC);

% Null for z-value
dPACnull = zeros(1,permutations);

% Permutate the signal
for permutation = 1:permutations

    % Show progress
    if mod(permutation,100) == 0
        fprintf('Combined sensors : Permutation %d / %d\n',permutation,permutations)
    end

    % -- cut-and-paste a random portion of the data; this preserves
    % -- temporal autocorrelation while removing the coupling
    cutLoc = 5 + randperm(length(phaseCombined)-10); % -- 5 and 10 prevent the first and last time points from being selected
    cutLoc = cutLoc(1);
    phaseShuffled = phaseCombined([cutLoc:end 1:cutLoc-1]);

    % Compute surrogate dPAC
    dPACnull(permutation) = abs(mean((exp(1i*phaseShuffled) - mean(exp(1i*phaseShuffled))) .* amplitudeCombined));
end

% dPAC z-score
meandPACnull = mean(dPACnull);
stddPACnull = std(dPACnull);
dPACz = (abs(dPAC) - meandPACnull) ./ stddPACnull;

% Store result
Zresults(NcombinedSensors + 1) = dPACz;

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
polarplot([0,angle(dPAC)],[0,abs(dPAC)],'k' , 'LineWidth',4)

saveas(figure(1),'PolarPlot','png')
hold off

%% Export results

% Write CSV textfiles with dPAC results (4 numbers, for sensors 1-3 and the
% combined one), one document for dPAC values, and one for z values.
csvwrite('000_GroupResults_dPAC.txt',dPACresults);
csvwrite('000_GroupResults_z.txt',Zresults);

% Stop timer
toc

% Notify user with sound
beep