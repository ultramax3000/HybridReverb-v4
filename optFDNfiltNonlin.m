function [gainsNonlin] = optFDNfiltNonlin(T60bands,delayTimes,centerFreqs,shelvingFreqs,...
    R,fs,numControlFreqs,initGain)
fftLen = 2^16;

numFreq = length(centerFreqs) + length(shelvingFreqs);
shelvingOmega = hertz2rad(shelvingFreqs, fs);
centerOmega = hertz2rad(centerFreqs, fs);

 
% Convert Frequency dependant T60 times to attenuation per sample (Schlecht
% EQ (2)-(5) & (27).

for i = 1:length(delayTimes)
    targetGains(:,i) = -60*(1./(fs*T60bands))*delayTimes(i,:);
end

% control frequencies are spaced logarithmically
controlFrequencies = logspace(log10(1), log10(fs/2.1),numControlFreqs+1);

% Interpolate for target magnitude response via attenuation gain per sample
 targetFreqs = [1, centerFreqs fs];
 
 for j=1:size(targetGains,2)
 tau(:,j) = interp1(targetFreqs, targetGains(:,j), controlFrequencies)';
 end

 
%% Desgin prototype of the filter magnitude responses
prototypeGain = 1; % dB
prototypeGainArray = prototypeGain * ones(numFreq+1,1);
prototypeSOS = proportionalParametricEQ(centerOmega, shelvingOmega, R, prototypeGainArray);
[B,prototypeH,prototypeW] = probeSOS (prototypeSOS, controlFrequencies, fftLen, fs);
B= B/prototypeGain;

% Plot the Magnitude response of the Filters
% figure; 
% semilogx(prototypeW,mag2db(abs(prototypeH)))
% ylim([0 11])
% xlim([10 fs/2])
% grid on;
% title('Prototype Magnitude Response (Interaction Matrix)')
% xlabel('Frequency [Hz]')
% ylabel('Magnitude [dB]')
 
%% Compute optimal parametric EQ gains (nonlinear solution)

% Initial gain in dB
%x0 = initGain*ones(numFreq+1, 1); 

% Bounds for the gains (between -10 and 10 dB)
upBound = 100*ones(numFreq+1, 1); 
lowBound = -upBound;

%optGainNonlin = lsqnonlin(nonlinFun,x0,upBound,lowBound);
options = optimoptions('fmincon','SpecifyObjectiveGradient',true); % indicate gradient is provided 

%Optimize gains for each delay-line
optGainNonlin = zeros(11,length(delayTimes));
for k = 1:size(tau,2)
    % Calculate x0 as the arithmetic mean of tau (32)
    x0 = zeros(numFreq+1, 1);
    x0(1) = mean(tau(:,k));
optGainNonlin(:,k) = fmincon(@(gamma) nonlinFun(gamma, B, tau(:,k)),x0,[],[],[],[],lowBound, upBound,[], options);
end

%finalCost = nonlinFun(optGainNonlin , B, T60bands,T60full');

%Plot approximated filter response
semilogx(controlFrequencies, B*optGainNonlin)

gainsNonlin = optGainNonlin;
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
 
figure
imagesc(B);
title(['Interaction Matrix with Filter Magnitude Responses (' num2str(numControlFreqs) ' Control Frequencies)']) 

end
