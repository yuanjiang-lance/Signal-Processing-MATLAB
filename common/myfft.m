function [f, ffty_Sig, fft_Sig] = myfft(Sig, SampFreq)
%
% Fast Fourier Transform of signals (single-side spectrum)
%
% ------- Input ---------
%  Sig: original signals (real or analytical), column listed
%  SampFreq: sampling frequency (Hz)
%
% ------- Output --------
%  f: Frequency bins (0: Fs/2)
%  ffty_Sig: the amplitude of FFT results, column listed
%  fft_Sig: FFT results, column listed
%
% Column listed means a group of signal or result lies in one column.
%
% Author: Yuan JIANG
% Time: 2023-07-30

NN = size(Sig, 1);
f = (0:NN-1)' * SampFreq/NN;
f = f(1: floor(NN/2));

if isreal(Sig)
    Sig = hilbert(Sig);
end

fft_Sig = fft(Sig)/NN;
fft_Sig = fft_Sig(1: floor(NN/2), :);
ffty_Sig = abs(fft_Sig);
