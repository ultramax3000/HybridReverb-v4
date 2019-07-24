function [T60full] = calcEDR100(signal,fs,frameSize,overlap,windowType,numControlFreqs)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%frameSize:   minimum frame length, in ms
%overlap:       fraction of frame overlapping
%windowType:    type of windowing used for each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Calculate 100 log-spaced frequencies between DC and Nyquist
controlFrequencies = logspace(log10(1), log10(fs/2.1),numControlFreqs+1);
freqs = controlFrequencies;

% calculate STFT frames
minFrameLen = fs*frameSize/1000;
frameLenPow = nextpow2(minFrameLen);
frameLen = 2^frameLenPow; % frame length = fft size
eval(['frameWindow = ' windowType '(frameLen);']);

[B,F,T] = spectrogram(signal,frameWindow,overlap*frameLen,freqs,fs);

[nBins,nFrames] = size(B);

B_energy = B.*conj(B);
B_EDR = zeros(nBins,nFrames);
for i=1:nBins
    B_EDR(i,:) = fliplr(cumsum(fliplr(B_energy(i,:))));
end
B_EDRdb = 10*log10(abs(B_EDR));

%  normalize EDR to 0 dB 
offset = max(max(B_EDRdb));
B_EDRdbN = B_EDRdb-offset;
minPlotDB = -60;

%Calculate the T60 times for each EDR Subband
for i=1:size(B_EDRdbN,1)
    T60full(:,i) = find(B_EDRdbN(i,:) <= -60,1, 'first');
end
T60full = T60full*(1-overlap)*frameLen/fs;

%% Plot
%truncate the plot below a given dB threshold
B_EDRdbN_trunc = B_EDRdbN;
for i=1:nFrames
    I = find(B_EDRdbN(:,i) < minPlotDB);
    if (I)
        B_EDRdbN_trunc(I,i) = minPlotDB;
    end
end
figure(gcf);clf;
surf(T,F/1000,B_EDRdbN_trunc,'EdgeAlpha',0.1);
set(gca, 'YScale','log', 'FontName','Palatino','FontSize',13);
view(130,30);
%title(['Normalized Energy Decay Relief (EDR) for ' num2str(numControlFreqs) ' log-spaced Bands']);
title('Normalized Energy Decay Relief (EDR)');
xlabel('Time (s)');ylabel('Frequency (kHz)');zlabel('Magnitude (dB)');
axis tight;zoom on;

end

