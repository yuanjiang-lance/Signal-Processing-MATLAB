

clear, clc, close all
path = '../../common';
addpath(path);
%% Signals
Fs = 1500;
T = 1;
t = 0: 1/Fs: T-1/Fs;

Sig1 = exp(1j*(2*pi*(300*t.^2)));
Sig2 = exp(1j*(340*pi*t + 150*sin(2*pi*t)));
Sig3 = exp(1j*(680*pi*t - 150*sin(2*pi*t)));
Sig4 = exp(1j*(2*pi*(600*t - 300*t.^2)));
Sig = Sig1 + Sig2 + Sig3 + Sig4;

sigshow(Sig, Fs);
[Spec, f] = mySTFT(Sig, Fs, 1024, 128);
stftshow(t, f, Spec);

%% Ridge detection
num = 4;
bw1 = Fs / 60;
delta = 20;
orderIF1 = 50;
lambda = 5;
orderIA1 = round(bw1*length(Sig)/Fs);

index = zeros(num, length(Sig));
reSig = Sig;
for i = 1: num
    Spec = mySTFT(reSig, Fs, 1024, 128);
    indextemp = DFRE(Spec, delta);
    [dSig, IFfit] = ICCD(Sig, Fs, f(indextemp), orderIA1, lambda, orderIF1);
    reSig = reSig - dSig;
    
end

figure
set(gcf,'Position',[20 100 640 500]);	    
set(gcf,'Color','w'); 
plot(t,f(index),'linewidth',3);
xlabel('Time / Sec','FontSize',24,'FontName','Times New Roman');
ylabel('Frequency / Hz','FontSize',24,'FontName','Times New Roman');
set(gca,'FontSize',24)
set(gca,'linewidth',2);
axis([0 T 0 Fs/2]);

%%
thrf = length(f) / 30;
index = RPRG(index, thrf);
figure
set(gcf,'Position',[20 100 640 500]);	    
set(gcf,'Color','w'); 
plot(t,f(findex),'linewidth',3);
xlabel('Time / Sec','FontSize',24,'FontName','Times New Roman');
ylabel('Frequency / Hz','FontSize',24,'FontName','Times New Roman');
set(gca,'FontSize',24)
set(gca,'linewidth',2);