function [f, ffty_Sig, fft_Sig] = myfft(Sig, SampFreq, figure)
% Calculating the FFT spectrum of a signal
% ------- Input ---------
%  Sig: Original Signal (real or analytical), column listed
%  SampFreq: Sampling Frequency
%  figure: 'on' = plot figure, otherwise do nothing
% ------- Output --------
%  f: Frequency bins (0: Fs/2)
%  ffty_Sig: FFT后的频谱幅值
%  fft_Sig: 信号的FFT结果

NN = size(Sig, 1);
f = (0:NN-1)' * SampFreq/NN;
f = f(1: floor(NN/2));

if isreal(Sig)
    Sig = hilbert(Sig);
end

fft_Sig = fft(Sig)/NN;
fft_Sig = fft_Sig(1: floor(NN/2), :);
ffty_Sig = abs(fft_Sig);

if nargin>2 && strcmp(figure, 'on')
    figure;
    set(gcf, 'position', [750 400 450 250], 'color', 'w');
    plot(f, ffty_Sig, 'b-', 'linewidth', 1);
    xlabel('Frequency/Hz', 'FontName', 'Times New Roman');
    ylabel('Amplitude/g', 'FontName', 'Times New Roman');
    set(gca, 'FontSize', 14, 'LineWidth', 1, 'FontName', 'Times New Roman');
end
