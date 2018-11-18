

function TruccoloData
%import Data
load MG29S1.mat;
L=length(LFP);
y=LFP;
x=t;

%Plot
figure
plot(x,y);
hold on;
title('Local Field Potential');
xlabel('Time (s)');
ylabel('Voltage (mV)');

%Plot Beginning and End of Seizure
plot([0,0], [-1000 600], 'r','LineWidth',1)
plot([66.5,66.5], [-1000 600], 'r','LineWidth',1)
xlim([-20,90]);
hold off

fft()


end