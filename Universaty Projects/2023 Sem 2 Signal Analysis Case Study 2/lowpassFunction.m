function lowpass = lowpassFunction(R,C)
%  LOWPASSFUNCTION Function takes the input of a capcitor and resistor, and creates a simple 
% first order low pass filter 
s=tf('s');
lowpass =1/(R*C*s+1);
end