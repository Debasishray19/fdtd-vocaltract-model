close all;
clc; clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEFINE TIME CONSTANTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SECOND = 1;    % Unit of time
MILLISECOND= 1e-3 * SECOND;
SRATE = 44100*15;                 % Sampling rate of sound card
dt    = 1/SRATE;               % Time period
total_audio_time = 1;   % Temporal length of sound generation
t = 0:dt:total_audio_time-dt;  % Temporal time steps
t_milli = t.*1000;
STEPS = length(t);

output_ug   = zeros(1, STEPS);
[airParam, vf_structuralParam, vf_flowParam, vf_matParam] = vf_SetVocalFoldParams();

for time_counter = 1:STEPS
    
    [vf_flowParam] = vf_TwoMassModel(SRATE, airParam, vf_structuralParam, vf_flowParam);
    
    output_ug(time_counter) = vf_flowParam.ug_next;   
end

% Plot the vocal fold volume velocity after normalization
norm_ug = output_ug./max(output_ug);
plot(t_milli(1:33075), norm_ug(1:33075), 'LineWidth', 1.5);
ylim([0 1.5]);
xlabel('Time [ms]');
ylabel('Normalized Volume Velocity [m^3/s]');
set(gcf, 'color', 'w');