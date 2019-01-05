Fs = 32000; % samples per second
dt = 1/Fs;
N = 1000; % Number of samples
data = csvread('/homes/nambigapathy/Desktop/TestData/DataFor_Matlab/tst.csv');
time = (0:1:(N-1))*dt;
timedata = data(:,2);
figure;
plot(time,timedata);
df = 1/(N*dt); % frequency increment
Nyq = 1/(dt*2); % Nyquist Frequency
freq = -Nyq:df:Nyq-df;
freqdata = zeros(size(timedata));
for i = 1 : N
for j = 1 : N
    freqdata(i) = freqdata(i) + timedata(j)*exp(-1i*2*pi*freq(i)*time(j));
end
end
mx = abs(freqdata);
figure;
plot(freq,real(freqdata),freq,imag(freqdata));