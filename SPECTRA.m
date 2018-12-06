%% Load Data
clc
clear all
%% Parameters 
N = length(velocitydata);       
T = 4.369;
dt = T/131071;
Fs = 1/dt;
time =  [0:4.369/131071:4.369];
F = [-N/2:(N-1)/2]./(N*dt);
meanU = mean(mean(velocitydata))
fluc = velocitydata- meanU;

%% Variance calculation 1
variance_direct = mean(var(fluc))

%% Calculate Spectra
S = zeros(size(fluc));
B = S;
for n =1:100
    X = fftshift(fft(fluc(:,n)));
    temp = (X.*conj(X))/(T)*dt^2;
    S(:,n) = temp;
    B(:,n) = autocorr(fluc(:,n),N-1);
end
%%
[Pxx , f] = pwelch(fluc,(2^17),0,[],Fs,'onesided','psd');
loglog(F,mean(S,2),'b')
hold on
loglog(f,mean(Pxx,2)/2,'r')
hold off

%% Variance 2 
ran = 10;
bd = [(N/2-ran+1):N/2,N/2+2:(N/2+ran)];
p = polyfit(F(bd)',mean(S(bd,:),2),2);
clf
plot(F,mean(S,2),'o')
hold on
plot(F(bd),polyval(p,F(bd)))
hold off

xlim([-10,10])
mS = mean(S,2);
mS(N/2+1) = p(end);
variace_spectra_1 = trapz(F,mS)
S0_1 = p(end)

%% Variance 3
[S1,F1] = pwelch(fluc,hamming(2^15),[],[],Fs);
ran =3;
bd = [2:ran];
p = polyfit([-flipud(F1(bd));F1(bd)],[flipud(mean(S1(bd,:),2));mean(S1(bd,:),2)],4);
clf
plot(F1,mean(S1,2))
hold on
plot([-flipud(F1(bd));0;F1(bd)],polyval(p,[-flipud(F1(bd));0;F1(bd)]),'-o')
hold off
xlim([-10,10])
variace_spectra_2 = trapz(F1,[p(end);mean(S1(2:end,:),2)])
S0_2 = p(end)/2
%%
integral_length_scale1 = (S0_1)/(2*variance_direct)
integral_length_scale2 = (S0_2)/(2*variance_direct)
integral_length_scale3 = max(cumsum(mean(B,2)).*dt)

