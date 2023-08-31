function sigshow(Sig,Fs)
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

figure
set(gcf,'Position',[20 100 640 500]);	    
set(gcf,'Color','w'); 
plot(t,Sig,'b','linewidth',2);
xlabel('Time (s)','FontSize',24,'FontName','Times New Roman');
ylabel('Amplitude (AU)','FontSize',24,'FontName','Times New Roman');
set(gca,'FontSize',24)
set(gca,'linewidth',2);
set(gca,'FontName','Times New Roman')