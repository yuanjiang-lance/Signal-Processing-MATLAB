function sigshow(Sig, Fs)
%
% Showing one-dimentional signal
%
% ---------- Input ----------
%  Sig: original signal
%  Fs: sampling frequency (Hz)
%
% Author: Yuan JIANG
% Time: 2023-08-31

N = length(Sig);
t = (0: N-1)/Fs;
Sig = real(Sig);

figure
set(gcf,'Position',[20 100 415 200]);	    
set(gcf,'Color','w'); 
plot(t,Sig,'b','linewidth',0.5);
xlabel('Time (s)','FontSize',12,'FontName','Times New Roman');
ylabel('Amplitude (AU)','FontSize',12,'FontName','Times New Roman');
set(gca,'FontSize',12)
set(gca,'linewidth',1);
set(gca,'FontName','Times New Roman')