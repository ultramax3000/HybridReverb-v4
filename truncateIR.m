function [truncTime, truncIR] = truncateIR(IR, fs, winSize)
%__________________________________________________________________________
%[truncation time, truncated array] =  truncate(array with IR, samplerate, window size in ms)

%Determines the end of the early reflections and truncates
%the IR. Output time of truncation in msec and truncated
%array are output.
%__________________________________________________________________________

%Function requires at least 3 input arguments
if nargin < 3
    error('Not enough input arguments.');
end
truncTime = []; 
truncIR = [];

%Window incoming IR file
window = fs * winSize/1000;
sizeIR = length(IR);

iterations = floor(sizeIR(1)/window);
time = linspace(0, sizeIR(1)/fs, iterations); %time increments

for i=1:iterations
    winIR = IR(window * (i-1)+ 1:window * i);
    
    %Find STANDARD DEVIATION and MEAN of window
    sigma = std(winIR);
    mu = mean(winIR);
    
    %Calculate percentage of samples outside 1 standard deviation
    outlier = 0;
    
    for n=1:window
        if winIR(n) < (mu - sigma) | winIR(n) > (mu + sigma)
            outlier = outlier + 1;
        end
    end
    
    %Store in array
    percent1(i) = outlier/window;
end

%Check for the first threshold crossing over 30% during the first 40 to
%60ms
met30 = [];
for i=8:-1:3
    if percent1(i) >= .3
        met30 = i;
    end
end

%Find local maximum beyond 60ms after first crossing of threshold
[~, truncIndex] = max(percent1(met30:met30+1));

truncTime = 1000*time(truncIndex + met30-1);
tempIR = IR(1:ceil(truncTime/1000 * fs)+32);

%Window the end of the IR with the second half of a hann window 
w = hann(64);
ERlength = length(tempIR);
tempIR((ERlength-31):ERlength) = (w(33:64))' .* tempIR((ERlength-31):ERlength);
%tempIR((ERlength-31):ERlength) = w(33:64) .* tempIR((ERlength-31):ERlength);

%Return truncated IR
truncIR = tempIR;

end

