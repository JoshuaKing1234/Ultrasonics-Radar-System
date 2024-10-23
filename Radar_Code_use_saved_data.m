%% Clear variables, close figures, and clear command window
clear all;
close all;
clc;

%% Define constants and parameters
SpeedSoundWave_ms = 343;            % [m/s]  -> Speed of sound
Fc_Hz = 40e3;                       % [Hz]   -> Carrier frequency
Fs = 192e3;                         % [Hz]   -> Sampling rate
Ts = 1/Fs;                          % Sampling period
lambda = SpeedSoundWave_ms/Fc_Hz;   % Wavelength


%% Walking
%filename = 'Saved_data\walk_towards_slow.mat';
%filename = 'Saved_data\walk_towards_fast';
%filename = 'Saved_data\walk_away_slow';
%filename = 'Saved_data\walk_away_fast';

%% Ball
%filename = 'Saved_data\ball_towards_slow';
%filename = 'Saved_data\ball_towards_medium';
filename = 'Saved_data\ball_towards_fast';

%% Car 
%filename = 'Saved_data\car_10';
%filename = 'Saved_data\car_20';
%filename = 'Saved_data\car_30';
%filename = 'Saved_data\car_40';


%%Other
%filename = "hand_back&forward.mat";
%filename = "hand_forward&back.mat";
%filename = "hand_forward_slow.mat";
%filename = "person_fast_medium.mat";
%filename = "person_forward_medium.mat";
%filename = "person_forward_slow.mat";
%filename= "Saved_data\max.mat";






load(filename); 
t = (0:length(RX_signal_rec)-1) / Fs;

%% Apply band-stop (notch) filter around 40 kHz
Order_bandstop = 4;                      % Filter order
Fc_lower = 39.8e3;                       % Lower cutoff frequency
Fc_higher = 40.2e3;                      % Higher cutoff frequency
[b_NF, a_NF] = butter(Order_bandstop, [Fc_lower/(Fs/2) Fc_higher/(Fs/2)], 'stop');
RX_signal_rec_NF = filtfilt(b_NF, a_NF, RX_signal_rec);

%% Baseband shifting
RX_signal_i = RX_signal_rec_NF .* sin(2*pi*Fc_Hz*t);
RX_signal_Q = RX_signal_rec_NF .* cos(2*pi*Fc_Hz*t);

%% Design and apply a low-pass filter
Fc_lowpass = 20e3;                     % [Hz] Cutoff frequency
Order = 10;                            % Filter order
[lp_b, lp_a] = butter(Order, Fc_lowpass/(Fs/2), 'low');  % Design low-pass filter
RX_signal_i_filtered = filtfilt(lp_b, lp_a, RX_signal_i);
RX_signal_Q_filtered = filtfilt(lp_b, lp_a, RX_signal_Q);

%% Combine the filtered signals into a complex signal
RX_signal = complex(RX_signal_i_filtered, RX_signal_Q_filtered);



% Spectrogram of the received signal
frame_length = 4096 ;%1024;                   % Frame length for spectrogram
%overlap = frame_length / 2;            % 50% overlap between frames
%overlap = round(frame_length * 0.6);  % 60% overlap
overlap = round(frame_length * 0.8);  % 80% overlap
%overlap = round(frame_length * 0.9);  % 90% overlap
%overlap = round(frame_length * 0.98);  % 90% overlap          % 50% overlap between frames
window = hamming(frame_length)';       % Hamming window
num_frames = floor((length(RX_signal) - frame_length) / (frame_length - overlap)) + 1;
spectrogram_matrix = zeros(frame_length, num_frames);
timevector_s = (0:length(RX_signal)-1) * (1/Fs);
TimeAxis_s = zeros(1, num_frames);

% Calculate FFT for each frame and fill the spectrogram matrix
for i = 1:num_frames
    start_index = (i - 1) * (frame_length - overlap) + 1;
    end_index = start_index + frame_length - 1;
    if end_index > length(RX_signal), break; end
    frame = RX_signal(start_index:end_index);          % Extract frame
    windowed_frame = frame .* window;                  % Apply window
    fft_frame = fftshift(fft(windowed_frame));         % Perform FFT and shift
    spectrogram_matrix(:, i) = fft_frame(1:frame_length);  % Store FFT result
    TimeAxis_s(i) = (timevector_s(start_index) + timevector_s(end_index))/2;
end

%% Frequency axis for spectrogram plot
FrequencyAxis_Hz = (-frame_length/2:1:(frame_length/2-1)) * (Fs/frame_length);

%% Plot spectrogram of received signal
figure;
imagesc(TimeAxis_s, FrequencyAxis_Hz, 20*log10(abs(spectrogram_matrix)));
xlabel('Time (s)', 'fontsize', 12);
ylabel('Frequency (Hz)', 'fontsize', 12);
title('Spectrogram of Received Signal', 'fontsize', 12);
axis xy;                               % y-axis increases from bottom to top
colormap('jet');
colorbar;

%% Plot the velocity spectrogram (zoomed-in)
maxSpeed_m_s = 10;
minSpeed_m_s = -10;
VelocityAxis_Hz = FrequencyAxis_Hz * lambda / 2; 
Velocity_Idx = find((VelocityAxis_Hz <= maxSpeed_m_s) & (VelocityAxis_Hz >= minSpeed_m_s));
VelocityAxis_Hz_subset = VelocityAxis_Hz(Velocity_Idx);
spectrogram_matrix_SubSet = spectrogram_matrix(Velocity_Idx, :);
spectrogram_matrix_SubSet = abs(spectrogram_matrix_SubSet) / max(max(abs(spectrogram_matrix_SubSet))); % Normalize

%% Plot zoomed-in velocity spectrogram
figure;
clims = [-50 0];
imagesc(TimeAxis_s, VelocityAxis_Hz_subset, 20*log10(abs(spectrogram_matrix_SubSet)), clims);
xlabel('Time (s)', 'fontsize', 12);
ylabel('Velocity (m/s)', 'fontsize', 12);
title('Spectrogram of Received Signal', 'fontsize', 12);
grid on;
colorbar;
colormap('jet');
axis xy;
