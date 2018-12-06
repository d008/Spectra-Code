%%
% This code has been graciously donated to you by Matt Fu. He wrote this
% code for his advanced turbulence course in 2016.
%% Load Data
clc
clear all
close all
load('velocitydata.mat');
%% Parameters 
N = length(velocitydata);       
T = 4.369; % period in seconds 
dt = T/(N-1); % time step in seconds 
Fs = 1/dt; % sampling frequency in Hz
time =  [0:T/N:T];
F = [-N/2:(N-1)/2]./(N*dt);
meanU = mean(mean(velocitydata));
fluc = velocitydata - meanU;

%% Variance Calculation 1
variance_direct = mean(var(fluc))

%% Calculate Spectra
S = zeros(size(fluc));
B = S;
for n =1:size(velocitydata,2)
    X = fftshift(fft(fluc(:,n)));
    temp = (X.*conj(X))/(T)*dt^2;
    S(:,n) = temp;
    B(:,n) = autocorr(fluc(:,n),N-1);
end
%% Plotting both Spectra and Computing Spectra from Pwelch
figure;
[Pxx , f] = pwelch(fluc,(2^17),0,[],Fs,'onesided','psd');
loglog(F,mean(S,2),'b')
hold on
loglog(f,mean(Pxx,2)/2,'r')
hold off
xlabel('Frequency [Hz]');
ylabel('PSD [m^2/s]');
legend('FFT','Pwelch');
grid on;
variance_spectra_1 = trapz(f, mean(Pxx,2));

%% Computing Variance from FFT
figure;
ran = 10;
bd = [(N/2-ran+1):N/2,N/2+2:(N/2+ran)];
p = polyfit(F(bd)',mean(S(bd,:),2),2);
clf
plot(F,mean(S,2),'o')
hold on
plot(F(bd),polyval(p,F(bd)))
hold off
grid on;
xlabel('Frequency [Hz]');
ylabel('FFT(u'') [m^2/s]');
legend('Data','Fit');

xlim([-10,10])
mS = mean(S,2);
mS(N/2+1) = p(end);
variance_spectra_2 = trapz(F,mS)
S0_1 = p(end);

%% Computing the Variance using another Pwelch
figure;
[S1,F1] = pwelch(fluc,hamming(2^15),[],[],Fs);
ran = 3;
bd = [2:ran];
p = polyfit([-flipud(F1(bd));F1(bd)],[flipud(mean(S1(bd,:),2));mean(S1(bd,:),2)],4);
clf
plot(F1,mean(S1,2))
hold on
plot([-flipud(F1(bd));0;F1(bd)],polyval(p,[-flipud(F1(bd));0;F1(bd)]),'-o')
hold off
xlim([-10,10])
grid on;
xlabel('Frequency [Hz]');
ylabel('FFT(u'') [m^2/s]');
legend('Data','Fit');
variance_spectra_3 = trapz(F1,[p(end);mean(S1(2:end,:),2)])
S0_2 = p(end)/2;
%% Computing the Integral Length Scales
integral_length_scale1 = (S0_1)/(2*variance_direct);
integral_length_scale2 = (S0_2)/(2*variance_direct);
integral_length_scale3 = max(cumsum(mean(B,2)).*dt);