%% Clear variables, close figures, and clear command window
clear all;
close all;
clc;

%% Define constants and parameters
SpeedSoundWave_ms = 343;            % [m/s]  -> Speed of sound
Fc_Hz = 40e3;                       % [Hz]   -> Carrier frequency
TimeDuration_s = 7;                 % [s]    -> Signal duration
Fs = 192e3;                         % [Hz]   -> Sampling rate
Ts = 1/Fs;                          % Sampling period
t = 0:Ts:(TimeDuration_s-Ts);       % Time vector for pulse
N = length(t);                      % Number of samples
lambda = SpeedSoundWave_ms/Fc_Hz;   % Wavelength

%% Generate the transmit signal
TxSignal = sin(2*pi*Fc_Hz*t);       % 40 kHz sine wave

%% Play out transmit signal through the speakers
soundsc(TxSignal, Fs, 24);          % Transmit the signal

%% Plot the transmit signal
figure; 
axes('fontsize', 12);
subplot(3,1,1);                     % Subplot 1 of 3
plot(t, TxSignal);                  % Plot transmit signal
xlabel('Time (s)', 'fontsize', 12);
xlim([0 0.001]);                    % Limit x-axis to first 1 ms
ylabel('Amplitude (linear)', 'fontsize', 12);
title('Transmit Signal', 'fontsize', 12);
grid on;

%% FFT and Frequency Domain Analysis
Y = fft(TxSignal);                  % FFT of the signal
Y_shifted = fftshift(Y);            % Shift zero-frequency component to center
f_fft_shifted = (-N/2:N/2-1)*(Fs/N); % Frequency vector with negative frequencies
Y_mag_shifted = abs(Y_shifted);     % Magnitude of shifted FFT
Y_dB_shifted = 20*log10(Y_mag_shifted + eps); % Magnitude in dB

%% Plot the signal in the frequency domain
subplot(3,1,2);                     % Subplot 2 of 3
plot(f_fft_shifted, Y_dB_shifted, 'LineWidth', 1.5);
xlabel('Frequency (Hz)', 'fontsize', 12);
ylabel('Magnitude (dB)', 'fontsize', 12);
title('Frequency Domain Spectrum', 'fontsize', 12);
grid on;
axis tight;

%% Record received samples from the microphone
RecLength_samples = length(TxSignal);
RecLength_s = RecLength_samples / Fs;
recObj = audiorecorder(Fs, 24, 1);
recordblocking(recObj, RecLength_s);    % Record audio for the same duration
RX_signal_rec = transpose(getaudiodata(recObj));  % Store recorded audio signal

% Save recorded data (uncomment to save)
 save_directory = 'Saved_data\max';
save(save_directory, 'RX_signal_rec');

%% Simulated data 
 %v = 10; % m/s
%Fd_Hz = 2*v/lambda;
%disp(Fd_Hz); 
 %RX_signal_sim = 1*sin(2*pi*Fc_Hz*t) + 1.5; % No hand
 %RX_signal_sim = sin(2*pi*Fc_Hz*t) + 1*sin(2*pi*(Fc_Hz + Fd_Hz)*t) + 1.5; % Hand moving
 %RX_signal_rec = RX_signal_sim;
 %RX_signal_i = RX_signal_rec .* sin(2*pi*Fc_Hz*t);
 %RX_signal_Q = RX_signal_rec .* cos(2*pi*Fc_Hz*t);



%% Walking
%filename = 'Saved_data\walk_towards_slow.mat';
%filename = 'Saved_data\walk_towards_fast';
%filename = 'Saved_data\walk_away_slow';
%filename = 'Saved_data\walk_away_fast';

%% Ball
%filename = 'Saved_data\ball_towards_slow';
%filename = 'Saved_data\ball_towards_medium';
%filename = 'Saved_data\ball_towards_fast';

%% Car 
%filename = 'Saved_data\car_10';
%filename = 'Saved_data\car_20';
%filename = 'Saved_data\car_30';
%filename = 'Saved_data\car_40';


%%Other
filename = "hand_back&forward.mat";
%filename = "hand_forward&back.mat";
%filename = "hand_forward_slow.mat";
%filename = "person_fast_medium.mat";
%filename = "person_forward_medium.mat";
%filename = "person_forward_slow.mat";
%load(filename); 
%t = (0:length(RX_signal_rec)-1) / Fs;

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
Fc_lowpass = 40e3;                     % [Hz] Cutoff frequency
Order = 10;                            % Filter order
[lp_b, lp_a] = butter(Order, Fc_lowpass/(Fs/2), 'low');  % Design low-pass filter
RX_signal_i_filtered = filtfilt(lp_b, lp_a, RX_signal_i);
RX_signal_Q_filtered = filtfilt(lp_b, lp_a, RX_signal_Q);

%% Combine the filtered signals into a complex signal
RX_signal = complex(RX_signal_i_filtered, RX_signal_Q_filtered);

%% Plot the received signal
subplot(3,1,3);                        % Subplot 3 of 3
plot(t, (RX_signal_rec));               % Plot the received signal
xlabel('Time (s)', 'fontsize', 12);
xlim([0 0.001]);                       % Limit x-axis to first 1 ms
ylabel('Amplitude (linear)', 'fontsize', 12);
title('Received Signal', 'fontsize', 12);
grid on;

%% Spectrogram of the received signal
frame_length = 4096 ;%1024;                   % Frame length for spectrogram
%overlap = frame_length / 2;            % 50% overlap between frames
%overlap = round(frame_length * 0.6);  % 60% overlap
overlap = round(frame_length * 0.8);  % 80% overlap
%overlap = round(frame_length * 0.9);  % 90% overlap
%overlap = round(frame_length * 0.98);  % 90% overlap

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
maxSpeed_m_s = 15;
minSpeed_m_s = -15;
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

