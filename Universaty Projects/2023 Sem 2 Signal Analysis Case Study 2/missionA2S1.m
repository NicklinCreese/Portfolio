%% EGB242 Assignment 2, Section 1 %%
% This file is a template for your MATLAB solution to Section 1.
%
% Before starting to write code, generate your data with the startHere
% script as described in the assignment task.

%% Initialise workspace
clear all; close all;
load DataA2 audioMultiplexNoisy fs sid;

% Begin writing your MATLAB solution below this line.

% Setup
samples = length(audioMultiplexNoisy);
T = samples / fs;

% Time vector
t = linspace(0, T, samples + 1); t(end) = [];

% Frequency vector
f = linspace(-fs/2, fs/2, samples + 1); f(end) = [];
f_kHz = f / 1000; % kHz

% Fourier transform
AudioMultiplexNoisy = fftshift(fft(audioMultiplexNoisy)) / fs;

% 1.1) Plot multiplexed audio signal in time and frequency domain
figure; sgtitle('Multiplexed Audio Signal');
subplot(2, 1, 1);
plotTime(t, audioMultiplexNoisy);
subplot(2, 1, 2);
plotMag(f_kHz, AudioMultiplexNoisy);

%% Demodulate and filter each carrier frequency
f_cutoff = 3e3;  % 3 kHz cut-off

% Locate carrier frequency peaks with mag > 5
[~, fc_peaks] = findpeaks(abs(AudioMultiplexNoisy), f, 'MinPeakHeight', 5);
fc_peaks = fc_peaks(fc_peaks >= 0); % Remove negative frequencies

% Demodulate and filter carrier frequencies
[filteredAudio, FilteredAudio] = demodulateAll(t, audioMultiplexNoisy, fs, f_cutoff, fc_peaks);

% Correct phase shift in 40.1 kHz signal
demodulatedAudio = audioMultiplexNoisy .* sin(2*pi*fc_peaks(3).*t);
filteredAudio(3, :) = lowpass(demodulatedAudio, f_cutoff, fs);
FilteredAudio(3, :) = fftshift(fft(demodulatedAudio)) / fs;

%% 1.2) Plot demodulated audio signals in time/frequency domain
for i = 1:length(fc_peaks)
    figure; sgtitle(['Demodulated Audio (f_c= ', num2str(fc_peaks(i)/1000), ' kHz)']);
	subplot(3, 1, 1);
	plotTime(t, filteredAudio(i, :));
	subplot(3, 1, 2);
	plotMag(f_kHz, FilteredAudio(i, :));
	subplot(3, 1, 3);
	plotPhase(f_kHz, FilteredAudio(i, :));
end

%% Play demodulated audio
n = 3; % Enter index of carrier frequency 1-5 [8120, 24020, 40100, 56290, 72080]
sound(filteredAudio(n, :), fs);

%% 1.3) Model frequency-dependent distortion
% Setup impulse
ts = 1 / fs; % Time per sample
impulse = [1 / ts, zeros(1, samples - 1)];

% Send impulse through channel
y = channel(sid, impulse, fs);

% Frequency response of channel
H = fftshift(fft(y)) / fs;

% Plot frequency response H(f) and magnitude spectrum of audio on same axes
figure; hold on;
plotMag(f_kHz, H);
plotMag(f_kHz, AudioMultiplexNoisy);
title('Frequency Response H(f) vs Magnitude Spectrum of AudioMultiplexNoisy', 'FontSize', 14);
legend('H(f)', 'AudioMultiplexNoisy');

%% 1.4) Reverse distortion and demodulate audio
% Reverse distortion
AudioMultiplexDenoised = AudioMultiplexNoisy ./ H;
audioMultiplexDenoised = ifft(ifftshift(AudioMultiplexDenoised)) * fs;

% Demodulate and filter denoised carrier frequencies
[filteredAudio2, FilteredAudio2] = demodulateAll(t, audioMultiplexDenoised, fs, f_cutoff, fc_peaks);

%% Plot de-noised multiplexed audio in time/frequency domain
figure; sgtitle('Denoised Multiplexed Audio Signal');
subplot(2, 1, 1);
plotTime(t, audioMultiplexDenoised);
subplot(2, 1, 2);
plotMag(f_kHz, AudioMultiplexDenoised);

%% Play demodulated/denoised audio
k = 5; % Enter index of carrier frequency 1-5 [8120, 24020, 40100, 56290, 72080]
sound(filteredAudio2(k, :), fs);

%% Prepare to Detone audio
% Identify min/max magnitude of tones on plot
figure; sgtitle('Denoised Multiplexed Audio Signal Before Detoning');
plotMag(f_kHz, AudioMultiplexDenoised);
xlim([-80 80]);
ylim([0 1]);
% Reduce vertical size of plot to fit title
plotPosition = get(gca, 'Position');
plotPosition(4) = 0.75;
set(gca, 'Position', plotPosition);

%% Detone, demodulate, filter, and normalise audio
% Set tone mag range
minToneMag = 0.3; 
maxToneMag = 0.7;

% Filter peaks within mag range
[mag, tone_peaks] = findpeaks(abs(AudioMultiplexDenoised), 'MinPeakHeight', minToneMag);
tone_peaks = tone_peaks(mag < maxToneMag);

% Detone audio
AudioMultiplexDenoisedDetoned = AudioMultiplexDenoised;
AudioMultiplexDenoisedDetoned(tone_peaks) = 0;
audioMultiplexDenoisedDetoned = ifft(ifftshift(AudioMultiplexDenoisedDetoned)) * fs;

% Demodulate and filter denoised/detoned carrier frequencies
[filteredAudio3, FilteredAudio3] = demodulateAll(t, audioMultiplexDenoisedDetoned, fs, f_cutoff, fc_peaks);

% Normalise audio
for i = 1:length(fc_peaks)
    filteredAudio3(i, :) = filteredAudio3(i, :) / max(abs(filteredAudio3(i, :)));
end

%% Plot denoised/detoned multiplexed audio in time/frequency domain
figure; sgtitle('Denoised and Detoned Multiplexed Audio Signal');
subplot(2, 1, 1);
plotTime(t, audioMultiplexDenoisedDetoned);
subplot(2, 1, 2);
plotMag(f_kHz, AudioMultiplexDenoisedDetoned);
xlim([-80 80]);
ylim([0 11]);

%% Play demodulated/denoised/detoned audio
j = 3; % Enter index of carrier frequency 1-5 [8120, 24020, 40100, 56290, 72080]
sound(filteredAudio3(j, :), fs);
