%% Execution Script for HybridReverb-v4

clear all; clc; close all

%% Specify files for IR and Dry Signal

audioInput = 'vox10.wav';
%[audio, fs] = audioread(audioInput);

%Use an Impulse:
audio = [0; 1; 0];
fs = 44100;

%Here, both a mono- and stereo IR have to be specified,because analysis and 
%optimization only make sense to be performed on a mono IR

%IR Mono:
%IRinputMono = 'academy_yard_mono.wav';     
IRMono = 'nave_cathedral_mono.wav';   
%IRinputMono = 'phipps_hall_huddersfield_mono.wav'; % EDR doesn't work (?)
%IRinputMono = 'BM7_Medium_Chamber_mono.wav';

IRStereo = 'nave_cathedral.wav';      %IR Stereo

%Read audio-files into arrays
[IR, fs,] = audioread(IRMono);
[IRstereo, fs2] = audioread(IRStereo);

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

% FDN Delay Times 16x16
%delayTimes = [89 97 107 113 149 211 263 293 401 421 433 443 577 601 641 661];
%delayTimes = [193 373 421 499 569 617 677 751 823 907 929 947 971 991 1019 1039];
%delayTimes = [443 1949 4409 5417 6421 7537 8863 9049 10799 11177 12791 13679 14891 15287 16339 17657];
%delayTimes = prime_power_delays(fs,16,1,100);

 delayTimes = round([10 11.6356 13.4567 16.734501 20.186199 25.741699... 
     31.469299 38.294399 46.6838 55.456699 65.175499 76.824303 88.562302 101.278... 
     115.397003 130.501999].*(fs/1000));

% FDN Delay Times 8x8
%delayTimes = round([1.39 2.51 11.1 19.4 41.2 60.2 73.6 96.6].*(fs/1000));
%delayTimes = prime_power_delays(fs,8,1,200);

%% Early Reflections Processing
%Truncate the IR to only the Early Reflections, using a 20ms Window:
[truncTime, truncIRL] = truncateIR(IRstereo(:,1), fs, 20);
[~, truncIRR] = truncateIR(IRstereo(:,2), fs, 20);

%CONVOLVE the IR with the dry audio:
dryWet = 100;
convERL = convolve(truncIRL, audio, dryWet);
convERR = convolve(truncIRR, audio, dryWet);

convER = [convERL,convERR];

%% RIR Anlysis
%Analyze EDR of IR to retrieve target subband T60 times:
frameSize = 30;         % minimum frame length, in ms
overlap = 0.75;         % fraction of frame overlapping (between 0-1)
windowType = 'hann';    % type of windowing used for each frame
freqs = [shelvingFreqs(1),centerFreqs,shelvingFreqs(2)];
numControlFreqs = 100;

[T60bands] = calcEDR(IR,fs,frameSize,overlap,windowType,numControlFreqs);

%Calculate global T60 acording to ISO 3382-1:2009:
[T60]=iosr.acoustics.irStats(IRMono,'graph',true,'spec','mean');

%% Optimization
%Optimize parametric filter gains from magnitude response (Linear Solution)
[gainsLin]   = optFDNfiltLin(T60bands,delayTimes',centerFreqs, shelvingFreqs,R,fs,numControlFreqs);

%Optimize parametric filter gains from T60 times (Nonlinear Solution)
%[gainsNonlin]   = optFDNfiltNonlin(T60full,delayTimes', centerFreqs, shelvingFreqs,R,fs,numControlFreqs);

%% Pass gain coefficients to FDN and process input signal

%16x16 FDN
FDN16multiOut = FDN16multi(audio,fs,centerFreqs,shelvingFreqs,R,gainsLin,delayTimes');

%8x8 FDN
%  FDNimp = FDN8(x,fs,centerFreqs,shelvingFreqs,R,gainsLin,delayTimes');

%% Mixing down 16x16 FDN to stereo
IACC = 0.6; %Interaural Correlation Coefficient between 0 and 1.0
numChan = 16;

FDNstereo = mixdown(FDN16multiOut,IACC,numChan);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Mix early and late reverberation to create final singal
%Window the late reflections
mix = 100;

FDNLength = length(FDNstereo);
startWind = (floor(truncTime/1000*fs)-32);
mixWindow = hann(64);
FDNstereo(1:startWind,:) = zeros(startWind,2);
FDNstereo((startWind+1):(startWind+32)) = FDNstereo((startWind+1):(startWind+32)).*mixWindow(1:32)';

%Zero-pad so that matrix dimensions agree
audioLength = length(audio);
convLength = length(convER);
audio((audioLength+1):FDNLength,1) = zeros((FDNLength - audioLength),1);
convER((convLength+1):FDNLength,1) = zeros((FDNLength - convLength),1);

%Add the wet audio before mixing with dry
wet(:,1) = (convER(:,1) + FDNstereo(:,1));
wet(:,2) = (convER(:,2) + FDNstereo(:,2));

%Add the early and late reflections and mix with dry audio and normalise
processedAudio(:,1) = (mix/100)*(wet(:,1)) + ((100-mix)/100)*audio;
processedAudio(:,2) = (mix/100)*(wet(:,2)) + ((100-mix)/100)*audio;

%% Demo against Convolution only.
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
% soundsc(processedAudio, fs);