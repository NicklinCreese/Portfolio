function highPassTransferFunction = highPassFunction(R, C)
%  HIGHPASSFUNCTION Create the transfer function of a simple RC high pass filter using component
% values R,C
s=tf('s');
highPassTransferFunction = R*C*s/(1+R*C*s);
end