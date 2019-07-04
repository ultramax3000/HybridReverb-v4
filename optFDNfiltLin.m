function [gains] = optFDNfiltLin(T60full,delayTimes,centerFreqs,shelvingFreqs,R,fs,numControlFreqs)
fftLen = 2^16;

numFreq = length(centerFreqs) + length(shelvingFreqs);
shelvingOmega = hertz2rad(shelvingFreqs, fs);
centerOmega = hertz2rad(centerFreqs, fs);

% control frequencies are spaced logarithmically
controlFrequencies = logspace(log10(1), log10(fs/2.1),numControlFreqs+1);

% target magnitude response via command gains
targetFreqs = [1, centerFreqs fs];
targetGains = [1; -1; 1; -1; 1; -1; 1; -1; 1; 1]*10; % dB
targetInterp = interp1(targetFreqs, targetGains, controlFrequencies)';

%% desgin prototype of the biquad sections
prototypeGain = 1; % dB
prototypeGainArray = prototypeGain * ones(numFreq+1,1);
prototypeSOS = proportionalParametricEQ(centerOmega, shelvingOmega, R, prototypeGainArray);
[B,prototypeH,prototypeW] = probeSOS (prototypeSOS, controlFrequencies, fftLen, fs);
% B = B / prototypeGain; % dB vs control frequencies ("Interaction Matrix")

% plot
figure(2);
semilogx(prototypeW,mag2db(abs(prototypeH)))
ylim([0 11])
xlim([10 fs/2])
grid on;
title('Prototype Magnitude Response (Interaction Matrix)')
xlabel('Frequency [Hz]')
ylabel('Magnitude [dB]')

%% Compute optimal parametric EQ gains (linear solution)
% Either you can use a unconstrained linear solver or introduce gain bounds
% at [-10dB,+10dB] with guaranteed no deviation from the self-similarity
% property. The plot shows the deviation between design curve and actual
% curve.
dBbound = 10;
x0 = dBbound*ones(numFreq+1, 1);
%x0 = ones(numFreq+1, 1);

% upperBound = [inf, dBbound * prototypeGain * ones(1,numFreq)];
% lowerBound = -upperBound;

%Limit upper bounds to 0 dB instead
upperBound = zeros(numFreq+1, 1);
lowerBound = -[inf, dBbound * prototypeGain * ones(1,numFreq)];

%Calculate tau (27) for each dealy-line:
for i = 1:length(delayTimes)
    tau(:,i) = -60*(1./(fs*T60full))*delayTimes(i,:);
end

%Optimize gains for each individual delay-line
optGainLin = zeros(11,length(delayTimes));
for k = 1:size(tau,2)
    optGainLin(:,k) = lsqlin(B, tau(:,k), [],[],[],[], lowerBound, upperBound,[]);
end

% optG = G\targetInterp; % unconstrained solution
%optimalSOS = proportionalParametricEQ( centerOmega, shelvingOmega, R, optGainLin );

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
heatmap(B);
title(['Interaction Matrix with Filter Magnitude Responses (' num2str(numControlFreqs) ' Control Frequencies)'])

end
