%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INITIALIZE MATLAB ENVIRONMENT AND SET PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear;
clc;

srate_mul = 15;
fprintf("Enter the sample rate multiplier= %d\n", srate_mul);

dur = 1;
fprintf("Enter the simulation duration= %f\n", dur);

simulationType = 3;
fprintf("Choose simulation type [1-Open Space 2-Regular Tube 3-Vowel Sound]= %d\n", simulationType);

vowel = 1;
fprintf("Select tube type [0: \\no vowel\\ 1:\\a\\ 2:\\i\\ 3:\\u\\ 4:\\e\\ 5:\\o\\ 6:\\I\\]= %d\n", vowel);

cross_sectionType = 2;
fprintf("Enter tube cross_section type [1: circular 2:elliptical 3:square]= %d\n", cross_sectionType);

junctionType = 1;
fprintf("Vocal tract junctiontype [1-centric 2-eccentric]= %d\n", junctionType)
    
sourceType = 2;
fprintf("Select source model type [1:\\Sine\\ 2:\\Gaussian\\ 3:\\Impulse\\ 4:\\Gaussian white noise\\ 5:\\VF model\\] = %d\n", sourceType);

pmlSwitch = 0;
fprintf("Switch on the PML layers [1:ON 0:OFF]= %d\n", pmlSwitch);

baffleSwitch = 0;
fprintf("Switch on circular baffle [1-Yes 0-No]= %d\n", baffleSwitch);

rad = 2;
fprintf("Select mouth termination condition [1: Open mouth 2:Dirichlet Boundary]= %d\n", rad);

boundaryInterpolation = 0;
fprintf("Switch on area interpolation [1:ON 0:OFF]= %d\n", boundaryInterpolation);

plotting = 1;
fprintf('Activate simulation plot [1-ON 0-OFF]= %d\n', plotting);

saveAudioData = 1;
fprintf('Save audio data [1-YES 0-NO]= %d\n', saveAudioData);

fprintf('\n--------------------------\nReady to go!\n--------------------------\n');
eval('talkingtube');

fprintf('\n--------------------------\nAll done!\nBye Bye!\n--------------------------\n\n');