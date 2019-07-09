function [gains] = optFDNfiltLin(T60bands,delayTimes,centerFreqs,shelvingFreqs,R,fs,numControlFreqs)
fftLen = 2^16;

numFreq = length(centerFreqs) + length(shelvingFreqs);
shelvingOmega = hertz2rad(shelvingFreqs, fs);
centerOmega = hertz2rad(centerFreqs, fs);

%Calculate attenuation per sample (Schlecht "Accurate Reverberation..."  Eq.2) 
for i = 1:length(delayTimes)
    targetGains(:,i) = -60*(1./(fs*T60bands))*delayTimes(i,:);
end

% control frequencies are spaced logarithmically
controlFrequencies = logspace(log10(1), log10(fs/2.1),numControlFreqs+1);

% Interpolate for target magnitude response via command gains
 targetFreqs = [1, centerFreqs fs];
 
 for j=1:size(targetGains,2)
 tau(:,j) = interp1(targetFreqs, targetGains(:,j), controlFrequencies)';
 end
%% desgin prototype of the biquad sections
prototypeGain = 1; % dB
prototypeGainArray = prototypeGain * ones(numFreq+1,1);
prototypeSOS = proportionalParametricEQ(centerOmega, shelvingOmega, R, prototypeGainArray);
[B,prototypeH,prototypeW] = probeSOS (prototypeSOS, controlFrequencies, fftLen, fs);
B= B/prototypeGain;

% plot
% figure(2);
% semilogx(prototypeW,mag2db(abs(prototypeH)))
% ylim([0 11])
% xlim([10 fs/2])
% grid on;
% title('Prototype Magnitude Response (Interaction Matrix)')
% xlabel('Frequency [Hz]')
% ylabel('Magnitude [dB]')

%% Compute optimal parametric EQ gains (linear solution)
% Either you can use a unconstrained linear solver or introduce gain bounds
% at [-10dB,+10dB] with guaranteed no deviation from the self-similarity
% property. The plot shows the deviation between design curve and actual
% curve.
 dBbound = 10;

%Set bounds between -10 and 10 dB
upperBound = [inf, dBbound * prototypeGain * ones(1,numFreq)];
lowerBound = -upperBound;

%Limit upper bounds to 0 dB instead
%upperBound = zeros(numFreq+1, 1);
%lowerBound = -[inf, dBbound * prototypeGain * ones(1,numFreq)];

%Optimize gains for each individual delay-line
optGainLin = zeros(11,length(delayTimes));
for k = 1:size(tau,2)
    optGainLin(:,k) = lsqlin(B, tau(:,k), [],[],[],[], lowerBound, upperBound,[]);
end

gains = optGainLin;
%% plot
% figure(3); hold on; grid on;
%
% [hOpt,wOpt] = freqz(optimalSOS,fftLen,fs);
% plot(controlFrequencies, targetInterp);
% plot(wOpt,mag2db(abs(hOpt)))
% plot(controlFrequencies, B*optGainLin)
%
% set(gca, 'xScale', 'log')
% ylim([-12 12])
% xlim([10 fs/2])
% title('Approximation Magnitude Response')
% legend('Target', 'Actual EQ', 'Design EQ','Location','SouthEast');

figure(4)
imagesc(B);
title(['Interaction Matrix with Filter Magnitude Responses (' num2str(numControlFreqs) ' Control Frequencies)'])

end
