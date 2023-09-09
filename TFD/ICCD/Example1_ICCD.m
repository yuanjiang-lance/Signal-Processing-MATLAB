

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
lambda1 = 5;
orderIA1 = 20;

index = zeros(num, length(Sig));
reSig = Sig;
for i = 1: num
    Spec = mySTFT(reSig, Fs, 1024, 128);
    indextemp = DFRE(Spec, delta);
    [dSig, IFfit] = ICCD(Sig, Fs, f(indextemp), orderIA1, lambda1, orderIF1);
    reSig = reSig - dSig;
    for j = 1: length(reSig)
        [~, indextemp(j)] = min(abs(f - IFfit(1,j)));
    end
    index(i, :) = indextemp;
end

figure
set(gcf,'Position',[20 100 415 350]);	    
set(gcf,'Color','w'); 
plot(t,f(index),'linewidth',2);
xlabel('Time (s)','FontSize',14,'FontName','Times New Roman');
ylabel('Frequency (Hz)','FontSize',14,'FontName','Times New Roman');
set(gca,'FontSize',14, 'FontName','Times New Roman', 'linewidth', 1);

%% Ridge path regrouping
thrf = 50;
index = RPRG(index, thrf);
figure
set(gcf,'Position',[20 100 415 350]);	    
set(gcf,'Color','w'); 
plot(t,f(index),'linewidth',2);
xlabel('Time (s)','FontSize',14,'FontName','Times New Roman');
ylabel('Frequency (Hz)','FontSize',14,'FontName','Times New Roman');
set(gca,'FontSize',14, 'FontName','Times New Roman', 'linewidth', 1);

%% Signal decomposition
lambda = 0.5;
orderIF = 5;
orderIA = 10;
[Sigest, IFest, IAest] = ICCD(Sig, Fs, f(index), orderIA, lambda, orderIF);

%% Final IFs
figure
set(gcf,'Position',[20 100 415 350]);	    
set(gcf,'Color','w'); 
plot(t,IFest,'linewidth',2);
xlabel('Time (s)','FontSize',14,'FontName','Times New Roman');
ylabel('Frequency (Hz)','FontSize',14,'FontName','Times New Roman');
set(gca,'FontSize',14, 'FontName','Times New Roman', 'linewidth', 1);

%% Reconstructed components
figure(501), hold on, box on
set(gcf,'Position',[20 100 415 200]);	    
set(gcf,'Color','w'); 
plot(t,IAest(1,:), 'r', 'linewidth',0.8);
plot(t,real(Sigest(1,:)), 'b', 'linewidth', 0.8);
plot(t,real(Sig1-Sigest(1,:)), 'k--', 'linewidth', 0.8);
xlabel('Time (s)'), ylabel('C1');
set(gca,'FontSize',12, 'FontName','Times New Roman', 'linewidth', 1);

figure(502), hold on, box on
set(gcf,'Position',[20 100 415 200]);	    
set(gcf,'Color','w'); 
plot(t,IAest(2,:), 'r', 'linewidth',0.8);
plot(t,real(Sigest(2,:)), 'b', 'linewidth', 0.8);
plot(t,real(Sig3-Sigest(2,:)), 'k--', 'linewidth', 0.8);
xlabel('Time (s)'), ylabel('C3');
set(gca,'FontSize',12, 'FontName','Times New Roman', 'linewidth', 1);

figure(503), hold on, box on
set(gcf,'Position',[20 100 415 200]);	    
set(gcf,'Color','w'); 
plot(t,IAest(3,:), 'r', 'linewidth',0.8);
plot(t,real(Sigest(3,:)), 'b', 'linewidth', 0.8);
plot(t,real(Sig2-Sigest(3,:)), 'k--', 'linewidth', 0.8);
xlabel('Time (s)'), ylabel('C2');
set(gca,'FontSize',12, 'FontName','Times New Roman', 'linewidth', 1);

figure(504), hold on, box on
set(gcf,'Position',[20 100 415 200]);	    
set(gcf,'Color','w'); 
plot(t,IAest(4,:), 'r', 'linewidth',0.8);
plot(t,real(Sigest(4,:)), 'b', 'linewidth', 0.8);
plot(t,real(Sig4-Sigest(4,:)), 'k--', 'linewidth', 0.8);
xlabel('Time (s)'), ylabel('C4');
set(gca,'FontSize',12, 'FontName','Times New Roman', 'linewidth', 1);