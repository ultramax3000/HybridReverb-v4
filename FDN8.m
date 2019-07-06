function [FDNout] = FDN8(audio,fs,centerFreqs,shelvingFreqs,...
    R,gains, delayTimes)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

x = audio;

segm = zeros(fs,1);
x = [x; segm; segm; segm];
dt = 1/fs;
tinput = 0:dt:(length(x)*dt)-dt;

y = zeros(1,length(x));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize FDN Matrix
% Hadamard Matrix
A = 1/sqrt(8) * hadamard(8);

% Householder Matrix (sometimes unstable)
% N = 8;
% u = ones(N,1);
% A = eye(N) -  (2 / N) * (u * u');

%Householder via Revursive Embedding (Experimental)
% matrix = 0.5*[1 -1 -1 -1; -1 1 -1 -1; -1 -1 1 -1; -1 -1 -1 1 ];
% A = 0.5*[matrix, -matrix; -matrix matrix];

% Noise Matrix
% N = 8;
% N = randn(N);
% [A,~] = qr(N);

% Diagonal Matrix
%A = diag(ones(1,8));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize Dely-Line buffers
m = delayTimes;

z1 = zeros(1,max(m));
z2 = zeros(1,max(m));
z3 = zeros(1,max(m));
z4 = zeros(1,max(m));
z5 = zeros(1,max(m));
z6 = zeros(1,max(m));
z7 = zeros(1,max(m));
z8 = zeros(1,max(m));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Proportional Parametric Filters
% Initialize parameters
Q = sqrt(R) / (R-1);

shelvingOmega   = hertz2rad(shelvingFreqs,fs);
centerOmega     = hertz2rad(centerFreqs,fs);

%% Calculate Coefficients for all filters
%Init buffer to store coefficients for each filter and each delay-line
EQcoefs = zeros(size(gains,1),6,length(delayTimes));

for j = 1:length(delayTimes)
    EQcoefs(:,:,j) = proportionalParametricEQ(centerOmega,shelvingOmega,R,gains(:,j));
end
% Plot the magnitude responses of the filters
% for l = 1:10
% freqz(EQcoefs(l,1:3),EQcoefs(l,4:6),fs,fs);
% hold on
% end

% Retain broadband gains and convert to linear
broadGains = 10.^(gains(1,:)/20);

%% FDN Loop
states = zeros(2,10,8);

for n = 1:length(y)
    temp = [z1(m(1)) z2(m(2)) z3(m(3)) z4(m(4)) z5(m(5)) z6(m(6)) z7(m(7)) z8(m(8))];
    
    y(n) = x(n) + z1(m(1)) + z2(m(2)) + z3(m(3)) + z4(m(4)) + z5(m(5)) + ...
        z6(m(6))+ z7(m(7)) + z8(m(8));
    
    %Apply Parametric Filters
    lastA = temp*A;
    
    %Call each delay-line and apply parametric filter
    for i = 1:length(lastA)
        % Loop through each filter within delay-line with individual coeff.  
        for k = 1:10
            [lastA(1,i), states(:,k,i)] = filter(EQcoefs(k,1:3,i),EQcoefs(k,4:6,i),lastA(1,i),states(:,k,i));
        end
    end
    
    %Apply broadband gains
    lastA = lastA.*broadGains;
    
    z1 = [(x(n) + lastA(1)) z1(1:length(z1)-1)];
    z2 = [(x(n) + lastA(2)) z2(1:length(z2)-1)];
    z3 = [(x(n) + lastA(3)) z3(1:length(z3)-1)];
    z4 = [(x(n) + lastA(4)) z4(1:length(z4)-1)];
    z5 = [(x(n) + lastA(5)) z5(1:length(z5)-1)];
    z6 = [(x(n) + lastA(6)) z6(1:length(z6)-1)];
    z7 = [(x(n) + lastA(7)) z7(1:length(z7)-1)];
    z8 = [(x(n) + lastA(8)) z8(1:length(z8)-1)];
    
end
FDNout = y;
%% Plot
dt = 1/fs;
plot(tinput,y,'g'); hold on;
plot(tinput,x,'k'); xlabel('Seconds'); ylabel('Amplitude');
title('8x8 Feedback Delay Network');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Spectrogram
figure(1)
subplot(2,1,1);
specgram(x);
title('Dry Signal');
subplot(2,1,2);
specgram(y);
title('Wet Signal');
end

