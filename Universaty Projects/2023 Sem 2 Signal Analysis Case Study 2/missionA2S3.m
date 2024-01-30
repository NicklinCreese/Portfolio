%%%% EGB242 Assignment 2, Section 3 %%
% This file is a template for your MATLAB solution to Section 3.
%
% Before starting to write code, generate your data with the startHere
% script as described in the assignment task.

%% Initialise workspace
clear all; close all;
load DataA2 imagesReceived;

% Begin writing your MATLAB solution below this line.
%% Supervisor Code // Display All Images // 
% Converting a received pixel stream to an image matrix
imageSize =[480,640];
imagesignalLength = imageSize(1,1)*imageSize(1,2);
imageNumber = 1;
im2D = reshape(imagesReceived(imageNumber,:),imageSize);
% Displaying an image in a figure
figure;
imshow(im2D);    
title(['Received Image of potential landing site ' ,num2str(imageNumber)]);
% Saving an image matrix as an image file (to include in report)
imwrite(im2D, 'filename.png');
%% 3.2 Visualise and identify image of potential landing site 1
fs= 1000; %Initilise fs 
siteToVis = 1
samples = length(imagesReceived);
T = samples/fs;
t = linspace(0,T,samples+1); t(end)=[];
f = linspace(-fs/2,fs/2,samples+1); f(end)=[];
IM1FFT=fftshift(fft(imagesReceived(siteToVis,:)))/fs;
figure()

subplot(3,1,1)
plot(t,imagesReceived(siteToVis,:))
title('Time Domain of image signal of landing site 1');
xlabel('Time [s]');
ylabel('Amplitude');
xlim([0,T]) 
ylim([2*min(rmoutliers(imagesReceived(siteToVis,:),"mean")),2*max(rmoutliers(imagesReceived(siteToVis,:),"mean"))]) % 2*min/2*max removing outliers 3 STD from mean

subplot(3,1,2)
plot(f,abs(IM1FFT))
title('Frequency Domain of image signal of landing site 1');
xlabel('Frequency [Hz]');
ylabel('Magnitude');
ylim([0,3.5]) % hardcoded limit to contrasts with later image

subplot(3,1,3)
plot(f,abs(IM1FFT))
title('Close up of frequency domain at point noise becomes significant');
xlabel('Frequency [Hz]');
ylabel('Magnitude');
xlim([110,170]) %hardcoded
ylim([0,1]) %hard coded
%% 3.3 TESTING ACTIVE FILTERS

R = 820; C =1*10^-6; % initalising R and C values
RC = R*C % Set R*C to RC for calculation simplification
s=tf('s') % initilse s variable


% Create transfer function for respective active filters and determine
% their characteristics
HActive1 = (1/RC^2)/(s^2+(2*s/RC)+1/(RC^2)) 
ltiH1=ltiview('bode',HActive1); %See affect of active filter, Is a lowpass, lowpass cutoff frequency of 125hz
obj = findobj(ltiH1, 'type', 'axes');
title(obj, 'Active Filter 1');

HActive2 = s^2/(s^2+(2*s/RC)+1/(RC^2))  
ltiH2=ltiview('bode',HActive2); %See affect of Filter, Is a highpass
obj = findobj(ltiH2, 'type', 'axes');
title(obj, 'Active Filter 2');


%% 3.3 High Pass and LowPass filter using R1,C1 and R2,C2 respectively
R1 = 1200;C1 = 10*10^-6; R2 = 1000; C2 = 4.7*10^-6; %Initilse component values

%Create Passive filter 1 transfer function using the conjunction of a
%highpass passive filter and low pass passive filter

HPassive1 = highPassFunction(R1,C1)*lowpassFunction(R2,C2); % create transfer function
ltiH1 = ltiview('bode',HPassive1) % Inspect characteristics of transfer function
obj = findobj(ltiH1, 'type', 'axes');
title(obj, 'Passive Filter 1');



% 3.3 Two Passive LowPass filters using R2,C2 and HighPass using R1,C1 

HPassive2 = lowpassFunction(R1,C1)*lowpassFunction(R2,C2); % create transfer function
ltiH2=ltiview('bode',HPassive2) % Inspect characteristics of transfer function, lowpass cutoff frequency of 11.7Hz
obj = findobj(ltiH2, 'type', 'axes');
title(obj, 'Passive Filter 2');

%% 3.4 Remove the noise from recieved signal and isolate 'clean' images
% imshow Selected Filter's (Active Filter 1) affect on image signal
yActive1=lsim(HActive1,imagesReceived(imageNumber,:),t); % simulate image 1 signal transfer through Passive filter 1 
image2dActive1 = reshape(yActive1,imageSize); % reshape it for display using imshow
figure()
imshow(image2dActive1)
title("Active Filter 1's affect on Image of Potential landing site 1")

%% 3.5 Display all clean landing site images and display image 1 Time and Frequency Domains of denoised image signal of landing site 1

FFTyActive1 = fftshift(fft(yActive1))/fs;
figure()
subplot(3,1,1)
plot(t,yActive1)
xlim([0,T]) 
ylim([2*min(rmoutliers(yActive1,"mean")),2*max(rmoutliers(yActive1,"mean"))]) % min/max removing outliers 3 STD from mean
title('Time Domain of denoised image signal of landing site 1');
xlabel('Time [s]');
ylabel('Amplitude');

subplot(3,1,2)
plot(f,abs(FFTyActive1))
title('Frequency Domain of denoised image signal of landing site 1');
xlabel('Frequency [Hz]');
ylabel('Magnitude');
ylim([0,3.5]) % hardcoded LIMIT to contrast

subplot(3,1,3)
plot(f,abs(FFTyActive1))
title('Close up of attenuation point (125Hz)');
xlabel('Frequency [Hz]');
ylabel('Magnitude');
xlim([110,170]) %attentuation point
ylim([0,1]) %attentuation point magnitude

for i=1:size(imagesReceived,1)
    im1DActive1=lsim(HActive1,imagesReceived(i,:),t);
    im2D = reshape(im1DActive1, imageSize);
    % Displaying an image in a figure
    figure;
    imshow(im2D);    
    title(['Recieved Denoised Image of potential landing site ' ,num2str(i)]);
    % Saving an image matrix as an image file (to include in report)
    imwrite(im2D, 'filename.png');
end
