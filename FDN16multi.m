function [FDNMultiOut] = FDN16multi(audio,fs,centerFreqs,shelvingFreqs,...
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
%A = 1/sqrt(16) * hadamard(16);

% Householder Matrix (Was found to give more "colourless" results)
 matrix = 0.5*[1 -1 -1 -1; -1 1 -1 -1; -1 -1 1 -1; -1 -1 -1 1 ];
 %Recursive Embedding
 A = 0.5*[matrix -matrix -matrix -matrix; -matrix matrix -matrix -matrix;...
    -matrix -matrix matrix -matrix; -matrix -matrix -matrix matrix];

% Noise Matrix
% N = 16;
% N = randn(N);
% [A,~] = qr(N);

% Diagonal Matrix
%A = diag(ones(1,16));
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

% Retain broadband gains and convert to linear
broadGains = 10.^(gains(11,:)/20);

%% FDN Loop
states = zeros(2,10,16);

for n = 1:length(y)
    temp = [z1(m(1)) z2(m(2)) z3(m(3)) z4(m(4)) z5(m(5)) z6(m(6)) z7(m(7)) z8(m(8))...
        z9(m(9)) z10(m(10)) z11(m(11)) z12(m(12)) z13(m(13)) z14(m(14)) z15(m(15)) z16(m(16))];
    
    y1(n) = x(n)+ z1(m(1));
    y2(n) = x(n)+ z2(m(2));
    y3(n) = x(n)+ z3(m(3));
    y4(n) = x(n)+ z4(m(4));
    y5(n) = x(n)+ z5(m(5));
    y6(n) = x(n)+ z6(m(6));
    y7(n) = x(n)+ z7(m(7));
    y8(n) = x(n)+ z8(m(8));
    y9(n) = x(n)+ z9(m(9));
    y10(n) = x(n)+ z10(m(10));
    y11(n) = x(n)+ z11(m(1));
    y12(n) = x(n)+ z12(m(1));
    y13(n) = x(n)+ z13(m(1));
    y14(n) = x(n)+ z14(m(1));
    y15(n) = x(n)+ z15(m(1));
    y16(n) = x(n)+ z16(m(1));
    
    
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
    z9 = [(x(n) + lastA(9)) z9(1:length(z9)-1)];
    z10 = [(x(n) + lastA(10)) z10(1:length(z10)-1)];
    z11 = [(x(n) + lastA(11)) z11(1:length(z11)-1)];
    z12 = [(x(n) + lastA(12)) z12(1:length(z12)-1)];
    z13 = [(x(n) + lastA(13)) z13(1:length(z13)-1)];
    z14 = [(x(n) + lastA(14)) z14(1:length(z14)-1)];
    z15 = [(x(n) + lastA(15)) z15(1:length(z15)-1)];
    z16 = [(x(n) + lastA(16)) z16(1:length(z16)-1)];
    
end
FDNMultiOut = [y1; y2; y3; y4; y5; y6; y7; y8; y9; y10; y11; y12; y13; y14; y15; y16];

end

