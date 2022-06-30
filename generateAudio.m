function  finalAudioOutput = generateAudio (Pr_Audio, srate, srate_mul)
    second = 1;
    hertz = 1/second;
    
    maxPressure = max(Pr_Audio);
    audioOutput = Pr_Audio./maxPressure;
    audioOutput(audioOutput>1) = 1;
    audioOutput(audioOutput<-1) = -1;
    
    % Generate a low pass butterwoth filter (IIR)
    fc = 22050*hertz; % Cutoff Frequency
    wc = (2*pi*fc)/srate;
    
    % create butterworth filter
    order = 2;  % filter order
    [b, a] = butter(order, wc);
    
    % Frequency response of the digital filter
    freqz(b, a);
    
    % filter the original signal
    filtered_audio = filter(b,a, audioOutput);
    
    % downsample the signal
    finalAudioOutput = filtered_audio(1:srate_mul:end);
    
    % play the audio
    sound(finalAudioOutput, 1*44100);
    
    % save audio file
    audiowrite('audioSound.mp4',finalAudioOutput,44100);
end