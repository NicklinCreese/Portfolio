%% EGB242 Assignment 2, Section 2 %%
% This file is a template for your MATLAB solution to Section 2.
%
% Before starting to write code, generate your data with the startHere
% script as described in the assignment task.
%% Initialise workspace
clear all; close all;
% Begin writing your MATLAB solution below this line.
% Express the transfer function symbolically.
s=tf('s');
% Formulate transfer function.
num = 1;
den = s*(s+0.5);
H1 = num/den;
% Generate time vector
samples = 10e4; T= 20;
t = linspace(0,T, samples + 1); t(end)=[];
% Produce an input step signal and record the resulting step response
inputStep = ones(size(t));
stepResponse = lsim(H1,inputStep,t);
figure;
subplot(2, 1, 1);
plot(t, inputStep);
xlabel('Time [s]');
ylabel('Input [V]')
title('Input step');
subplot(2, 1, 2);
plot(t, stepResponse);
xlabel('Time [s]');
ylabel('Output Yaw Angle [radians]');
title('Step Response');
%% Feedback circuit
% Create transfer function
Kp =1;
num = H1;
den = 1+ H1*Kp;
H2 = num/den;
% Create step response
stepResponse = lsim(H2,inputStep,t);
figure;
plot(t, stepResponse);
xlabel('Time [s]');
ylabel('Output Yaw Angle [radians]');
title('Step response');
%% Gain Feedack circuit changing Kfb
Kfwd = 1;
num = Kfwd*H1;
for Kfb=[0.1, 0.2, 0.5, 1, 2]
den = 1+ Kfwd*Kfb*H1*Kp;
H3 = num/den;
stepResponse = lsim(H3,inputStep,t);
plot(t, stepResponse);
hold on
title('Step response at the different Kfb values');
xlabel('Time [s]');
ylabel('Output Yaw Angle [radians]');
end
legend('0.1','0.2','0.5','1','2')
%% Gain Feedack circuit changing Kfwd
Kfb = 1;
for Kfwd=[0.1, 0.2, 0.5, 1, 2]
num = Kfwd*H1;
den = 1+ Kfwd*Kfb*H1*Kp;
H3 = num/den;
stepResponse = lsim(H3,inputStep,t);
plot(t, stepResponse);
hold on
xlabel('Time [s]');
ylabel('Output yaw angle [radians]');
title('Step response at the different Kfwd values');
end
legend('0.1','0.2','0.5','1','2')
%% Gain feedback circuit chosen Kfwd and Kfb
Kfb = 1/(2*pi);
Kfwd = 0.7596;
num = Kfwd*H1;
den = 1+ Kfwd*Kfb*H1*Kp;
cameraTF = num/den;
stepResponse = lsim(cameraTF,inputStep,t);
figure;
plot(t, stepResponse);
xlabel('Time [s]');
title('Step response at the chosen Kfb and Kfwd values');
%% Start Simulation
[startIm, finalIm] = cameraPan(1/12, 7/12, cameraTF);
