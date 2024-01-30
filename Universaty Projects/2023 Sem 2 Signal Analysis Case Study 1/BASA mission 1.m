%% EGB242 Assignment 1 %%
% This file is a template for your MATLAB solution.
%
% Before starting to write code, record your test audio with the record
% function as described in the assignment task.
%% Load recorded test audio into workspace
clear all; close all;
load DataA1;
% Begin writing your MATLAB solution below this line.
% Begining Noisy Signal
sound(audio,fs)
lenAudio =length(audio)
t1=linspace(0,10,lenAudio);
plot(t1,audio);
title('Audio amplitude verus Time')
ylabel('Audio amplitude')
xlabel('Time(s)')
%% Calculating And Using complex fourier series using hand calculated function
t=linspace(0,10,lenAudio);
orderOftransform = 5;
cn = zeros([2*orderOftransform+1,2]);
for n=0:orderOftransform
for j=1:2
n=n*((rem(j,2)*2)-1);
p1=i*pi()*n;
a1=((5*exp(-p1)-5*exp(-4))/(4-p1));
b1=(2*p1+3*exp(-p1)-5*p1*exp(-p1)-3)/((pi()*n)^2);
final1 = (a1+b1)/2;
if n==0
final1 = 2.3635;
end
cn(n+orderOftransform+1,1) = final1;
cn(n+orderOftransform+1,2) = p1;
end
end
sumSignal= (cn(:,1)).*exp(cn(:,2).*t);
nApprox = real(sum(sumSignal));
audioClean = audio-nApprox;
%% Comparison from Noisy Audio to Approximated Noise Signal
figure()
plot(t,audio)
hold on
plot(t,nApprox,LineWidth=1,LineStyle="--")
ylabel('Amplitude')
xlabel('Time(s)')
title('Noisy audio output compared to the approximated noise signal')
%% Plotting / listening to Cleaned Audio Signal
figure()
plot(t1,audioClean)
ylabel('Amplitude')
xlabel('Time(s)')
title('Microphone Input signal accounting for Noise')
sound(audioClean,fs)
%%
f = linspace(-fs/2,fs/2,lenAudio+1);
f(end) = [];
fc = 38500;
X = fftshift(fft(audioClean))/fs;
%%{
figure()
plot(f,abs(X));
title('Frequency magnitude spectrum of Audio signal'), xlabel('Frequency
(Hz)'),ylabel('Magnitude')
channelQuiet = channel(11042338,zeros(size(t)));
C = fftshift(fft(channelQuiet))/fs;
figure()
plot(f,abs(C));
title('Frequency magnitude spectrum of Quiet Channel'), xlabel('Frequency
(Hz)'),ylabel('Magnitude')
%}
%%{
mod = audioClean.*(cos(2*pi*fc*t));
M= fftshift(fft(mod))/fs;
simSignal = channel(11042338,mod);
%Plotting Frequency magnitude spectrum of modulated audio and close ups at
%microphone frequency
figure()
tiledlayout (2,2,'TileSpacing','loose')
nexttile([1 2])
S=fftshift(fft(simSignal))/fs;
plot(f,abs(S));
title({'Frequency magnitude spectrum of modulated audio',' signal passed through
simulated channel'}), xlabel('Frequency (Hz)'),ylabel('Magnitude')
nexttile(4);
plot(f,abs(S));
xlim([30000, 47000] )
ylim([0 0.025])
title('Close up on output signal'), xlabel('Frequency (Hz)'),ylabel('Magnitude');
nexttile(3);
plot(f,abs(S));
title('Close up on input signal'), xlabel('Frequency (Hz)'),ylabel('Magnitude')
xlim([-47000,-30000] )
ylim([0 0.025])
%}
% Demodulate Audio and apply lowpass filter to reconstruct audio signal
demodAudio= mod.*(cos(2*pi*fc*t));
T= fftshift(fft(demodAudio))/fs;
figure()
plot(f,abs(T));
title('Demodulated Signal after transmission'), xlabel('Frequency
(Hz)'),ylabel('Magnitude')
audioRecieved=lowpass(demodAudio,fc*2/fs);
sound(audioRecieved,fs)
T= fftshift(fft(audioRecieved))/fs;
figure()
plot(f,abs(T));
title('Audio recieved signal after low pass filter'), xlabel('Frequency
(Hz)'),ylabel('Magnitude')
%% Listening to Active Frequencies in Channel
ChosenSig=1 %1-3
SignalList=[19200 57600 76800]
demodQuiet = channelQuiet.*(cos(2*pi*SignalList(ChosenSig)*t));
T= fftshift(fft(demodQuiet))/fs;
figure()
plot(f,abs(T));
title('Demodulated Signal of Chosen Sig'), xlabel('Frequency
(Hz)'),ylabel('Magnitude')
audioRecieved1=lowpass(demodQuiet,SignalList(ChosenSig)*2/fs);
sound(audioRecieved1,fs);
%% Quantisation
fs2= 22050;
audioResampled = resample(audioRecieved,fs2,fs);
sound(audioResampled,fs2);
