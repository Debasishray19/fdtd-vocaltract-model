function [bwStore]  = calculateFormantPosBW(audioSignal, sourceSignal, srate)

    %Sample Rate
    fs = srate;
    
    % Flip the source signal so that both audioSignal and SourceSignal can
    % have size=1*M
    if size(audioSignal,1) ~= size(sourceSignal,1)
        sourceSignal=sourceSignal';
    end
    
    audiosignalfft = fft(audioSignal);
    sourceSignalfft = fft(sourceSignal);
    len_audioSignal = length(audioSignal);
    
    TransferFunction = audiosignalfft./sourceSignalfft(1:len_audioSignal);
    finalTransferFunction = 10*log10(abs(TransferFunction).^2); %yVal
        
    n = length(audioSignal)-1;  
    df = fs/n;
    f=0:df:fs;  % xVal
    
    numBandwidth = input('Number of Formants Bandwidth: ');
    bwStore = zeros(3, numBandwidth);
    bwCounter = 1;
    valCounter = 2;
    
    flag=0;
    while bwCounter <= numBandwidth
        while finalTransferFunction(valCounter) > finalTransferFunction(valCounter-1)...
                && finalTransferFunction(valCounter) > finalTransferFunction(valCounter+1)
            
            % Store the formants' coordinates
            xFormant = f(valCounter);
            yFormant = finalTransferFunction(valCounter);
            yBWFormant = yFormant - 3;
            bwStore(2, bwCounter) = xFormant;
            bwStore(3, bwCounter) = yFormant;
            % store xVals and yVals in an array [left Side of the array]
            formantArrayCreateCounterLeft = 1;
            while finalTransferFunction(valCounter - formantArrayCreateCounterLeft) > yBWFormant
                if finalTransferFunction(valCounter - formantArrayCreateCounterLeft) <...
                       finalTransferFunction(valCounter - formantArrayCreateCounterLeft-1)
                   break;
                end
                formantArrayCreateCounterLeft=formantArrayCreateCounterLeft+1;
            end
            xFormantLeft = f(valCounter - formantArrayCreateCounterLeft: valCounter);
            yFormantLeft = finalTransferFunction(valCounter - formantArrayCreateCounterLeft: valCounter);
            
            formantArrayCreateCounterRight = 1;
            while finalTransferFunction(valCounter + formantArrayCreateCounterRight) > yBWFormant
                formantArrayCreateCounterRight=formantArrayCreateCounterRight+1;
            end
            xFormantRight = f(valCounter: valCounter + formantArrayCreateCounterRight);
            yFormantRight = finalTransferFunction(valCounter: valCounter + formantArrayCreateCounterRight);
            
            xBWFormantLeft = interp1(yFormantLeft, xFormantLeft, yBWFormant, 'linear', 'extrap');
            xBWFormantRight = interp1(yFormantRight, xFormantRight, yBWFormant, 'linear', 'extrap');
            
    
            %Find BW
            BW = xBWFormantRight-xBWFormantLeft;
            
            %Store BW
            bwStore(1,bwCounter)=BW;
            
            % Go to the next one
            bwCounter=bwCounter+1;
            
            valCounter = valCounter+1;
            
            % switch on th flag to avoid increasing valCounter again
            flag = 1;
        end
        
        if flag == 0
            valCounter = valCounter+1;
        else
            flag=0;
        end
    end 
    
%     save('data_025.mat','f','finalTransferFunction', 'bwStore');
end