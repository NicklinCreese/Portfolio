function [filteredAudio, FilteredAudio] = demodulateAll(t, audio, fs, f_cutoff, fc_peaks)
  	% Create matrices
	demodulatedAudio = zeros(length(fc_peaks), length(audio));
	filteredAudio = zeros(length(fc_peaks), length(audio));

	for i = 1:length(fc_peaks)
		% Demodulate audio at each carrier frequency (fc)
		demodulatedAudio(i, :) = audio .* cos(2*pi*fc_peaks(i).*t);
		% Filter demodulated audio
		filteredAudio(i, :) = lowpass(demodulatedAudio(i, :), f_cutoff, fs);
	end

	% Fourier Transform entire matrix
	FilteredAudio = fftshift(fft(filteredAudio, [], 2)) / fs;
end
