function [processedSignal] = convolve(IR, audio, mix)
%__________________________________________________________________________
%[processedAudio] = pilot(IR array, audio array, mix wet%)
%
%Takes an array of an IR and an array of dry audio then convolves,
%returning the mixed array. Only supports a mono dry audio with
%mono IR.
%__________________________________________________________________________

% Function requires at least 3 input arguments
if nargin < 3
    error('Not enough input arguments!');
end

processedSignal = [];

%Check the number of channels in dry audio-file
[audioLength, audioChan] = size(audio);
%Check the number of channels on IR-file
[IRlength, IRChan] = size(IR);

if (audioChan > 1)
    error('Too many channels in dry audio, please choose a mono file!');
end

if (IRChan > 1)
    error('Too many channels in IR, please choose a mono file!');
end
    %Convolve IR with dry audio
    if(audioChan == IRChan)
        tempAudio = conv(audio, IR);
    end
    
    %Allocate zero-padded buffer to same size as convolved audio
    audio((audioLength + 1):(length(tempAudio))) = zeros((length(tempAudio)...
        -(audioLength)),1);
    
    %Return convolved mixed with dry audio to an array
    processedSignal = (mix/100)*tempAudio + (1-mix/100)*audio;
    
end

