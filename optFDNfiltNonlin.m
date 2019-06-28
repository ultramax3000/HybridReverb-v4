function [gainsNonlin] = optFDNfiltNonlin(T60full,delayTimes,centerFreqs,shelvingFreqs,...
    R,fs,numControlFreqs)
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
 
%% Compute optimal parametric EQ gains (nonlinear solution)
x0 = ones(numFreq+1, 1); % initial gain in dB
upBound = 10*ones(numFreq+1, 1); % bounds for the gains
lowBound = -upBound;

%optGainNonlin = lsqnonlin(nonlinFun,x0,upBound,lowBound);
options = optimoptions('fmincon','SpecifyObjectiveGradient',true); % indicate gradient is provided 

% Convert Frequency dependant T60 times to attenuation per sample (Schlecht
% EQ (2)-(5) & (27), taking the average of all delay-line times (experimental solution.

%Calculate tau (27) for each dealy-line: 
tau = zeros(length(T60full),length(delayTimes)); 
for i = 1:length(delayTimes)
    tau(:,i) = -60*(1./(fs*T60full))*delayTimes(i,:);
end

%Optimize gains for each delay.line
optGainNonlin = zeros(11,length(delayTimes));
for k = 1:size(tau,2)
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
 
figure(4)
heatmap(B);
title(['Interaction Matrix with Filter Magnitude Responses (' num2str(numControlFreqs) ' Control Frequencies)']) 

end
