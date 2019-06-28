function [twoChAudio] = mixdown(multiChAudio, IACC, numChan)
%__________________________________________________________________________
%[stereo output audio] = mixdown(array holding audio, IACC value, number of
% input channels)
%
%Mixes 16 channels to 2 uncorrelated channels, then mixes them according
%to the IACC value (Interaural Corellation). Returns a stereo array of audio.
%
%Rebecca Stewart
%August 2006
%__________________________________________________________________________

if numChan == 16
    mixArray = [1 -1 1 -1 1 -1 1 -1 1 -1 1 -1 1 -1 1 -1; -1 1 1 -1 -1 1 1 ...
        -1 -1 1 1 -1 -1 1 1 -1];
end

if numChan == 4
    mixArray = [1 -1 1 -1; -1 1 1 -1];
end

mixedAudio = mixArray * multiChAudio;

theta = asin(IACC)/2;

twoChAudio(:,1) = cos(theta) * mixedAudio(1,:) + sin(theta) * mixedAudio(2,:);
twoChAudio(:,2) = cos(theta) * mixedAudio(2,:) + sin(theta) * mixedAudio(1,:);

end