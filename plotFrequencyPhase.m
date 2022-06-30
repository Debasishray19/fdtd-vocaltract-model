function plotFrequencyPhase(audioSignal, sourceSignal, srate)
    
    % Define sampling rate
    fs = srate;
    
    % Flip the sourceSignal to match its size to the audioSignal
    if size(audioSignal,1) ~= size(sourceSignal,1)
        sourceSignal=sourceSignal';
    end
    
    audiosignalfft = fft(audioSignal);
    sourceSignalfft = fft(sourceSignal);
    len_audioSignal = length(audioSignal);
    
    TransferFunction = audiosignalfft./sourceSignalfft(1:len_audioSignal);
    finalTransferFunction = 10*log10(abs(TransferFunction).^2);
        
    audiosignalfft_phase = angle(audiosignalfft);
    n = length(audioSignal)-1;  
    df = fs/n;
    f=0:df:fs;
    
    f_eccen_cir = f;
    finalTransferFunction_eccen_cir = finalTransferFunction;
    save('circular_eccentric_a.mat','f_eccen_cir', 'finalTransferFunction_eccen_cir');
    
    figure;
    plot(f,finalTransferFunction, 'LineWidth', 1); % plot Fourier Transform
    title('Amplitude Spectrum Analysis');
    xlabel('Frequency [Hz]');
    ylabel('Amplitude [dB]');
    axis 'auto y';
    yLim = ylim();
    axis([2 10050 yLim(1) yLim(2)])
    
    figure
    plot(f,audiosignalfft_phase); % plot Fourier Transform
    title('Phase Spectrum Analysis');
    xlabel('Frequency [Hz]');
    ylabel('Phase [Radian]');
    axis 'auto y';
    yLim = ylim();
    axis([2 22050 yLim(1) yLim(2)])
end