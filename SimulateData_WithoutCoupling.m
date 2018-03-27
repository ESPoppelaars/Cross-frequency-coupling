clear
clc

% Number of files to create
nFiles = 20;

% Cycle through the files
for file = 1:nFiles
    
    % Show progress
    fprintf('File %d / %d\n',file,nFiles)

    %% Settings

    % Sensors to combine
    combinedSensors = [5 , 38 , 40];

    % Number of sensors
    nSensors = 64;

    % Number of epochs
    epochs = 6;

    % Duration of epoch
    epochTime = 8;

    % Sampling frequency
    fSample = 512;

    % Matrix with all signals
    signalMatrix = zeros(epochs * epochTime * fSample,nSensors);

    % Cycle through the sensors
    for sensor = 1:nSensors

        % If the sensor is in the sensor list
        if sum(sensor == combinedSensors)

            % Initialise signal variable
            signal = zeros(epochs * epochTime * fSample,1);

            % Cycle through the epochs
            for epoch = 1:epochs

                %% Delta signal

                % Gaussian variable
                g = .05;

                % Signal frequency
                %fDelta = 2.5;
                fDelta = rand*3 + 1;

                % End time of signal
                tEnd = epochTime;

                % Time series
                t = 1/fSample:1/fSample:tEnd;

                % Number of sample
                N = length(t);

                % Delta signal
                
                % Comment out to simulate data with bias:
                % Simulate data without bias:
                delta = cos(2*pi*fDelta * t);
                              
                %Create amplitude variability with a random frequency between 1 and 50 Hz.
                deltaAmplitude = 1 + sin((rand*49 + 1)*t * 2*pi + rand*2)/(2 + rand);
                
                % Multiply the generated signal with the amplitude.
                delta = delta.*deltaAmplitude;

                %% Beta Signal

                % Beta signal frequency
                %mean of fBeta = 22;
                fBeta = rand*16 + 14;

                % Beta signal
                
                %Simulate data without coupling:
                beta = cos(2*pi*fBeta * t);
                %Create amplitude variability with a random frequency between 1 and 50 Hz.
                betaAmplitude = 1 + sin((rand*49 + 1)*t * 2*pi + rand*2)/(2 + rand);
                % Multiply the generated signal with the amplitude  (only without coupling).
                beta = beta.*betaAmplitude;

                %% Noise

                % Noise amplitude
                %Medium noise:
                k = 6;
                %Large noise:
                %k = 20;

                % define the length of the vector
                % ensure that the M is even
                if rem(N,2)
                    M = N+1;
                else
                    M = N;
                end

                % White noise
                whiteNoise = k*randn(1, M);

                % Fast Fourier
                fourier = fft(whiteNoise);

                % prepare a vector for 1/f multiplication
                NumUniquePts = M/2 + 1;
                n = 1:NumUniquePts;
                n = sqrt(n);

                % multiply the left half of the spectrum so the power spectral density
                % is proportional to the frequency by factor 1/f, i.e. the
                % amplitudes are proportional to 1/sqrt(f)
                fourier(1:NumUniquePts) = fourier(1:NumUniquePts)./n;

                % prepare a right half of the spectrum - a copy of the left one,
                % except the DC component and Nyquist frequency - they are unique
                fourier(NumUniquePts+1:M) = real(fourier(M/2:-1:2)) -1i*imag(fourier(M/2:-1:2));

                % Inverse Fourier
                pinkNoise = ifft(fourier);

                % prepare output vector y
                pinkNoise = real(pinkNoise(1, 1:N));

                % ensure zero mean value
                pinkNoise = pinkNoise - mean(pinkNoise);

                %% Signal

                % Epoch signal
                subSignal = delta + beta + pinkNoise;

                % Append epoch signal to participant signal
                signal((epoch-1)*epochTime*fSample+1:epoch*epochTime*fSample) = subSignal';

                % Append signal to signal matrix
                signalMatrix(:,sensor) = signal;
            end
        end
    end

    %% Export data

   dlmwrite(strcat('SimulatedData',num2str(file),'.csv'),signalMatrix,',',1,0)

end

%Plot the last epoch of the last 'participant'.
plot(delta)
hold on
plot(beta)
hold off