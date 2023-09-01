function stftshow(t, f, Spec)
%
% Showing time-frequency representation (TFR), such as STFT
%
% ---------- Input ------------
%  t: time bins
%  f: frequency bins
%  Spec: time-frequency spectrum
%  
% Author: Yuan JIANG
% Time: 2023-08-31

figure
imagesc(t, f, abs(Spec));
JET = colormap(jet);
mycolor = cat(1, [1 1 1], JET(2:end,:));
colormap (JET);   
box on;
colorbar off;
set(gcf,'position',[846.6,340.2,414.4,364]);
set(gca,'linewidth',1.5,'fontsize',14,'fontname','Times New Roman');
xlabel('Time (s)','FontSize',14);
ylabel('Frequency (Hz)','FontSize',14);
axis([0 t(end)+0.005 0 f(end)+1]);  % constrain the TFR showing range, could be modified
view(0,90);
set(gca,'YDir','normal');
set(gcf,'Color','w');

end

