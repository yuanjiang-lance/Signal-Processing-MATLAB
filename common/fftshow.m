function fftshow(f, ffty_Sig)
%
% Showing frequency spectrum
%
% --------- Input -----------
%  f: frequency bins
%  ffty_Sig: FFT result or its amplitude
%
% Author: Yuan JIANG
% Time: 2023-07-30

if ~isreal(ffty_Sig)
    ffty_Sig = abs(ffty_Sig);
end

figure;
set(gcf, 'position', [750 400 450 250], 'color', 'w');
plot(f, ffty_Sig, 'b-', 'linewidth', 1);
xlabel('Frequency/Hz', 'FontName', 'Times New Roman');
ylabel('Amplitude/g', 'FontName', 'Times New Roman');
set(gca, 'FontSize', 14, 'LineWidth', 1, 'FontName', 'Times New Roman');