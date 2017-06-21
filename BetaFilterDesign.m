function BetaFilter = BetaFilterDesign
%BetaFilterDesign Returns a discrete-time filter object.

% MATLAB Code
% Generated by MATLAB(R) 8.5 and the Signal Processing Toolbox 7.0.
% Generated on: 24-Mar-2017 12:51:13

% Butterworth Bandpass filter designed using FDESIGN.BANDPASS.

% All frequency values are in Hz.
Fs = 128;             % Sampling Frequency

Fstop1 = 12;          % First Stopband Frequency
Fpass1 = 14;          % First Passband Frequency
Fpass2 = 30;          % Second Passband Frequency
Fstop2 = 32;          % Second Stopband Frequency
Astop1 = 24;          % First Stopband Attenuation (dB)
Apass  = 1;           % Passband Ripple (dB)
Astop2 = 24;          % Second Stopband Attenuation (dB)
match  = 'passband';  % Band to match exactly

% Construct an FDESIGN object and call its BUTTER method.
h  = fdesign.bandpass(Fstop1, Fpass1, Fpass2, Fstop2, Astop1, Apass, ...
                      Astop2, Fs);
BetaFilter = design(h, 'butter', 'MatchExactly', match);

% [EOF]
