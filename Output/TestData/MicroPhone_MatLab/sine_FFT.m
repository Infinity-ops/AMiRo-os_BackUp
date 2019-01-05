close all; clear all; clc;
Fs = 32000; % sampling frequency
data = csvread('/homes/nambigapathy/Desktop/TestData/DataFor_Matlab/tst.csv');
t = data(:,1);  % data
signal =  % Signal
nfft = 1024; %Length of FFT
x = Discrete_Fourier_Transform(signal);
x = x(1:nfft/2);
mx = abs(x);
f= (0:nfft/2-1)* Fs/nfft;

% Plot time-domain signal
subplot(2,1,1);
plot(t, signal);
ylabel('Amplitude'); xlabel('Time (secs)');
title(' Input Signal');


%TEST_Narmada
subplot(2,1,2);
plot(f,mx);
ylabel('Amplitude'); xlabel('Time (secs)');

title('Noisy Input Signal');


