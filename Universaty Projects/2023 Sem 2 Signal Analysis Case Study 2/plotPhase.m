function plotPhase(f, Audio)
	% plot(f, unwrap(angle(Audio)));
	plot(f, angle(Audio));
	title('Phase Spectrum');
	xlabel('Frequency [kHz]');
	ylabel('Phase [rad]');
end