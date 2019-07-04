%% Execution Script for HybridReverb-v3

clear all; clc; close all

%Specify files for IR and Dry Signal
audioInput = 'vox10.wav';

%IRinputMono = 'academy_yard_mono.wav';     %IR Mono
IRinputMono = 'nave_cathedral_mono.wav';   %IR Mono
%IRinputMono = 'phipps_hall_huddersfield_mono.wav'; % EDR doesn't work (?)
%IRinputMono = 'BM7_Medium_Chamber_mono.wav';

%IRinputStereo = 'nave_cathedral.wav';      %IR Stereo

%Read audio-files into arrays
[IR, fs,] = audioread(IRinputMono);
[audio, fs2] = audioread(audioInput);

%Make sure that sampling frequencies are equal for IR and dry signal
if (fs > fs2 || fs < fs2)
    error('IR and input signal need to have same fs!');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Global Parameters
centerFreqs = [62.5, 125, 250, 500, 1000, 2000, 4000, 8000];
shelvingFreqs = [31.25, 16000];
R = 3.6;                %relates to the Q-Factor
Q = sqrt(R) / (R-1);

% FDN Delay Times
%delayTimes = [89 97 107 113 149 211 263 293 401 421 433 443 577 601 641 661];
%delayTimes = [193 373 421 499 569 617 677 751 823 907 929 947 971 991 1019 1039];
%delayTimes = [443 1949 4409 5417 6421 7537 8863 9049 10799 11177 12791 13679 14891 15287 16339 17657];

delayTimes = prime_power_delays(fs,16,1,100);

%% Early Reflections Processing
%Truncate the IR to only the Early Reflections, using a 20ms Window:
[truncTime, truncIR] = truncateIR(IR, fs, 20);

%CONVOLVE the IR with the dry audio:
dryWet = 100;
convER = convolve(truncIR, audio, dryWet);

%% RIR Anlysis
%Analyze EDR of IR to retrieve target subband T60 times:
frameSize = 30;         % minimum frame length, in ms
overlap = 0.75;         % fraction of frame overlapping (between 0-1)
windowType = 'hann';    % type of windowing used for each frame
freqs = [shelvingFreqs(1),centerFreqs,shelvingFreqs(2)];
numControlFreqs = 100;

[T60full] = calcEDR(IR,fs,frameSize,overlap,windowType,numControlFreqs);

%Calculate global T60 acording to ISO 3382-1:2009:
[T60]=iosr.acoustics.irStats(IRinputMono,'graph',true,'spec','mean');

%Calculate Delay-line Lengths according to global T60:
%delayTimes = DelayLineLengths(16,fs,T60);

%% Optimization
%Optimize parametric filter gains from magnitude response (Linear Solution)
[gainsLin]   = optFDNfiltLin(T60full,delayTimes',centerFreqs, shelvingFreqs,R,fs,numControlFreqs);

%Optimize parametric filter gains from T60 times (Nonlinear Solution)
%[gains]   = optFDNfilt(T60full,delayTimes, centerFreqs, shelvingFreqs,R,fs,numControlFreqs);
%[gainsNonlin]   = optFDNfiltNonlin(T60full,delayTimes, centerFreqs, shelvingFreqs,R,fs,numControlFreqs);

%% Pass gain coefficients to FDN and process input signal

%Excite FDN with an Impulse:
x = [0; 1; 0];
FDNimp = FDN16(x,fs,centerFreqs,shelvingFreqs,R,gainsLin,delayTimes');

%16x16 FDN:
%FDNout = FDN16(audio,fs,centerFreqs,shelvingFreqs,R,gainsLin,delayTimes);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Prepare sample for listening
% [inputSignal] = audioread(audioInput);
% 
% %Generate only convolved signal
% dryWetMix = 50;
% [inputIR, fs] = audioread(IRinputStereo);
% [convOnlyL] = convolve(inputIR(:,1), inputSignal, dryWetMix);
% [convOnlyR] = convolve(inputIR(:,2), inputSignal, dryWetMix);
% convOnly = [convOnlyL, convOnlyR];
% 
% soundsc(inputSignal, fs);
% pause(3.5)
% soundsc(convOnly, fs);
% pause(3.5)
% soundsc(FDNout, fs);