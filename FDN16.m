function [FDNout] = FDN16(audio,fs,centerFreqs,shelvingFreqs,...
    R,gains, delayTimes)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

x = audio;

segm = zeros(fs,1);
x = [segm; x; segm; segm];
dt = 1/fs;
tinput = 0:dt:(length(x)*dt)-dt;

%%
y = zeros(1,length(x));
b = 0.9*ones(1,16);
%b = rand(1,16);
c = 0.9*ones(1,16);
c = rand(1,16);
% Matrix Gain coefficient |g|<1
g = 0.239999999999999999999;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Puckette and Stautner Matrix
% ??????????
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% using Hadamard Matrix
A = g*(1/2)*hadamard(16);
%% using Householder Matrix
%matrix = 0.5*[1 -1 -1 -1; -1 1 -1 -1; -1 -1 1 -1; -1 -1 -1 1 ];
%A = 0.5*[matrix -matrix -matrix -matrix; -matrix matrix -matrix -matrix;...
%   -matrix -matrix matrix -matrix; -matrix -matrix -matrix matrix];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 16 Delay lines, use prime
m = delayTimes;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Delay lines
z1 = zeros(1,max(m));
z2 = zeros(1,max(m));
z3 = zeros(1,max(m));
z4 = zeros(1,max(m));
z5 = zeros(1,max(m));
z6 = zeros(1,max(m));
z7 = zeros(1,max(m));
z8 = zeros(1,max(m));
z9 = zeros(1,max(m));
z10 = zeros(1,max(m));
z11 = zeros(1,max(m));
z12 = zeros(1,max(m));
z13 = zeros(1,max(m));
z14 = zeros(1,max(m));
z15 = zeros(1,max(m));
z16 = zeros(1,max(m));
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

%% FDN Loop
states = zeros(2,10,16);

for n = length(segm):length(y)
    temp = [z1(m(1)) z2(m(2)) z3(m(3)) z4(m(4)) z5(m(5)) z6(m(6)) z7(m(7)) z8(m(8))...
        z9(m(9)) z10(m(10)) z11(m(11)) z12(m(12)) z13(m(13)) z14(m(14)) z15(m(15)) z16(m(16))];
    
    y(n) = x(n) + c(1)*z1(m(1)) + c(2)*z2(m(2)) ...
        + c(3)*z3(m(3)) + c(4)*z4(m(4)) + c(5)*z5(m(5)) + c(6)*z6(m(6)) ...
        + c(7)*z7(m(7)) + c(8)*z8(m(8)) + c(9)*z9(m(9)) + c(10)*z10(m(10)) ...
        + c(11)*z11(m(11)) + c(12)*z12(m(12)) + c(13)*z13(m(13)) + c(14)*z14(m(14)) ...
        + c(15)*z15(m(15)) + c(16)*z16(m(16));
    
    %Apply Parametric Filters
    lastA = temp*A;
    
    %Call each delay-line and apply parametric filter
    for i = 1:length(lastA)
        % Loop through each filter within that delay-line with individual
        % coefficients  
        for k = 1:10
            [lastA(1,i), states(:,k,i)] = filter(EQcoefs(k,1:3,i),EQcoefs(k,4:6,i),lastA(1,i),states(:,k,i));
        end
    end
    
    z1 = [(x(n)*b(1) + lastA(1)) z1(1:length(z1)-1)];
    z2 = [(x(n)*b(2) + lastA(2)) z2(1:length(z2)-1)];
    z3 = [(x(n)*b(3) + lastA(3)) z3(1:length(z3)-1)];
    z4 = [(x(n)*b(4) + lastA(4)) z4(1:length(z4)-1)];
    z5 = [(x(n)*b(5) + lastA(5)) z5(1:length(z5)-1)];
    z6 = [(x(n)*b(6) + lastA(6)) z6(1:length(z6)-1)];
    z7 = [(x(n)*b(7) + lastA(7)) z7(1:length(z7)-1)];
    z8 = [(x(n)*b(8) + lastA(8)) z8(1:length(z8)-1)];
    z9 = [(x(n)*b(9) + lastA(9)) z9(1:length(z9)-1)];
    z10 = [(x(n)*b(10) + lastA(10)) z10(1:length(z10)-1)];
    z11 = [(x(n)*b(11) + lastA(11)) z11(1:length(z11)-1)];
    z12 = [(x(n)*b(12) + lastA(12)) z12(1:length(z12)-1)];
    z13 = [(x(n)*b(13) + lastA(13)) z13(1:length(z13)-1)];
    z14 = [(x(n)*b(14) + lastA(14)) z14(1:length(z14)-1)];
    z15 = [(x(n)*b(15) + lastA(15)) z15(1:length(z15)-1)];
    z16 = [(x(n)*b(16) + lastA(16)) z16(1:length(z16)-1)];
    
end
FDNout = y;
%% Plot
dt = 1/fs;
plot(tinput,y,'g'); hold on;
plot(tinput,x,'k'); xlabel('Seconds'); ylabel('Amplitude');
title('16x16 Feedback Delay Network');
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

