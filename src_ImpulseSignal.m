
%% The src_ImpulseResponse will return an impulse exciation [band-pass filter]
%% Given Inputs:
% sample_rate   : audio sampling rate
% lcf           : low cut-off frequency
% hcf           : high cut-off frequency
% maxExcitation : excitation magnitude
% excitationX   : excitation velocity

function excitationX = src_ImpulseSignal(sample_rate, maxExcitation, lcf, hcf)
    
    excitationX = zeros(sample_rate, 1);
    % Define your cut-off frequency as fraction of sampling rate
    lcf = lcf/sample_rate;
    hcf = hcf/sample_rate;
    
    % roll-off: Steepness of transmission frequency between the pass band and stop band
    bw = 0.001;    
    M  = 4/bw; % Number of points
    plotFigure = 0;
    
    imp = zeros(M+1,1);
    lp_sinc = zeros(M+1,1);
    hp_sinc = zeros(M+1,1);
    
    % STEP1: Design a low pass filter using sinc function
    sumVal = 0;
    for i=1:M+1
        if i-1~=M/2
            lp_sinc(i) = sin(2*pi*hcf * ((i-1)-M/2)) / ((i-1)-M/2);
        else
            lp_sinc(i) = 2*pi*hcf;
        end
        % multiply by Hamming window
        lp_sinc(i) = lp_sinc(i) * (0.54 - 0.46*cos(2*pi*(i-1)/M)); 
        sumVal = sumVal+lp_sinc(i);
    end
       
    % STEP2: Normalize the value of low pass filter
    for i=1:M+1
        lp_sinc(i) = lp_sinc(i)/sumVal;
    end
    
    % STEP3: Design a high pass filter using sinc function
    % [Note]: We can design a high pass filter using low pass filter and 
    % spectral inversion. 
    % First designed the low pass filter using cut-off frequency.
    % Then invert the low pass filter and add 1 to the sample at the shape
    % of symmetry. (Chapter 14: The scientist and engineer's guide to signal 
    % processing)
    
    sumVal = 0; % to normilize sincs
    for i=1:M+1
        if i-1~=M/2
            hp_sinc(i) = sin(2*pi*lcf * ((i-1)-M/2)) / ((i-1)-M/2);
        else
            hp_sinc(i) = 2*pi*lcf;
        end
        hp_sinc(i) = hp_sinc(i) * (0.54 - 0.46*cos(2*pi*(i-1)/M)); % multiply by Hamming window
        sumVal = sumVal+hp_sinc(i);
    end
    
    % normlize and invert, to turn into high pass sinc
    for i=1:M+1
        hp_sinc(i) = hp_sinc(i)/sumVal;
        hp_sinc(i) = -hp_sinc(i);
    end
    hp_sinc(M/2+1) = hp_sinc(M/2+1) + 1;
    
    % STEP 4: Construct a window sync(band pass) filter
    for i=1:M+1
        imp(i) = hp_sinc(i) + lp_sinc(i);
        imp(i) = -imp(i);
    end
    imp(M/2+1) = imp(M/2+1) + 1;
        
    ex = imp*maxExcitation;
    
    if M+1 >= sample_rate
        excitationX  = ex;
    else
        excitationX(1:M+1) = ex;
    end 
    
    % STEP 5: Plot the figure
    if plotFigure==1
        
    end
end