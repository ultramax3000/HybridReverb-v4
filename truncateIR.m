function [truncTime, truncIR] = truncateIR(IR, Fs, winSize)
%__________________________________________________________________________
%[truncation time, truncated array] =
% truncate(array containing IR, sampling frequency, window size in ms)
%%
%Determines the end of the early reflections and truncates
%the IR. The output time of truncation in msec and truncated
%array are output.
%
%Rebecca Stewart
%July 2006
%__________________________________________________________________________
%Function requires at least 3 inputs argument
if nargin < 3
    error('Not enough input arguments.');
end
truncTime = []; truncIR = [];
%Window wav
window = Fs * winSize/1000;
sizeInput = length(IR);
iterations = floor(sizeInput(1)/window);
time = linspace(0,sizeInput(1)/Fs,iterations); %time increments

for i=1:iterations
    winIn = IR(window * (i-1)+1:window*i);
    %Find standard deviation and mean of window
    sigma = std(winIn);
    mu = mean(winIn);
    %Calculate percentage of samples outside 1 standard deviation
    out1 = 0;
    for n=1:window
        if winIn(n) < (mu - sigma) | winIn(n) > (mu + sigma);
            out1 = out1 + 1;
        end
    end
    %Store in array
    percent1(i) = out1/window;
end
%Check for first time over 30% in first 40 to 160 ms
for i=8:-1:3
    if percent1(i) >= .3
        met30 = i;
    end
end
%Find local max over the 60 ms after the first crossing of the threshold
[maxValue, truncIndex] = max(percent1(met30:met30+1));
truncTime = 1000*time(truncIndex+met30-1);
tempIR = IR(1:ceil(truncTime/1000 * Fs)+32);
%Hann window
w = hann(64);
%Window the end of the IR with the falling edge of window function
l = length(tempIR);
tempIR((l-31):l) = w(33:64) .* tempIR((l-31):l);
%Return truncated IR
truncIR = tempIR;
end