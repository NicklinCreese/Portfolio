function plotMag(f, Audio)
    plot(f, abs(Audio));
    title('Magnitude Spectrum');
    xlabel('Frequency [kHz]');
    ylabel('Magnitude');
end