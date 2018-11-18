function truccolodata
%Sumaiya Sayeed

close all; clear all;

% Import Data
load('MG29S1');
L=length(LFP);
y=LFP;
x=t;

%Plot
figure
plot(x,y); %plot EEG signal
hold on
title('Local Field Potential');
xlabel('Time [s]');
ylabel('LFP [mV]');
%plot beginning and end of seizure
plot([0,0],[-1000 600],'r','LineWidth',1)
plot([66.5,66.5],[-1000 600],'r','LineWidth',1)
xlim([-20 90])
hold off

% Normalized Spectrogram - uses built-in Matlab function
s = spectrogram(LFP);
figure('Name','Spectrogram');
spectrogram(LFP,'yaxis');

% FFT on entire signal
figure
Y=fft(LFP);
Fs = 2/(t(2)-t(1));   % Sampling frequency                    
T = 1/Fs;             % Sampling period       
sze = size(LFP); 
L=sze(2); % Length of signal 
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1); %P1 = y-axis of FFT
f = Fs*(0:(L/2))/L;
plot(f,P1)
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')

% Lowpass for FFT
lowsize=10000 %plots t
figure
plot(f(1:lowsize),P1(1:lowsize)) %plots only the first 10,000 data points (< 180 Hz)
title('Single-Sided Amplitude Spectrum of X(t) (lowpass)')
xlabel('f (Hz)')
ylabel('|P1(f)|')

% STFT to develop a spectrogram using a periodogram
dt = x(2)-x(1)
Fs = 1/dt
T = 0.5; %time window
df = 1/T; %frequency resolution
fNQ = Fs/2;

F = [0:df:fNQ] %frequency domain
n = Fs*T %number of samples in each time window
j = 0;
round(n/2)

for k=round(n+1):round(n/2):round(L-n)
    j=j+1;
    time(j)=x(k); 
    signal=LFP(k-n+1:k);
    signal=signal(:); %column vector
    signal=signal-mean(LFP); %normalizing
    signalf=fft(signal); 
    temp=dt^2*1/T*abs(signalf).^2; %temporary variable
    Spectrum(:,j)=10*log10(2*temp(1:n/2+1));
    
end

figure
imagesc(time,F,Spectrum);
title('Spectrogram Based on Periodogram');
xlabel('Time [s]')
ylabel('Frequency')
set(gca,'YDir','normal')
colorbar;
hold on
dim = [.3 .2 .2 .3]; 
str = 'Noice the yellow areas between time 40-60';
annotation('textbox',dim,'String',str,'FitBoxToText','on');

hold off

end