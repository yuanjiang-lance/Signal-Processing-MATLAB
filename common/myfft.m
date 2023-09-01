function [f, ffty_Sig, fft_Sig] = myfft(Sig, Fs)
%
% Fast Fourier Transform of signals (single-side spectrum)
%
% ------- Input ---------
%  Sig: original signals (real or analytical), column listed for
%       multi-signal, or a row/column vector for just one signal
%  Fs: sampling frequency (Hz)
%
% ------- Output --------
%  f: frequency bins (0: Fs/2)
%  ffty_Sig: the amplitude of FFT results, column listed
%  fft_Sig: FFT results, column listed
%
% Column listed means a group of signal or result lies in one column.
%
% Author: Yuan JIANG
% Time: 2023-07-30

if size(Sig, 1) == 1
    Sig = Sig.'; 
    rowVec = 1; 
else
    rowVec = 0;
end

NN = size(Sig, 1);
f = (0:NN-1)' * Fs/NN;
f = f(1: floor(NN/2));

if isreal(Sig)
    Sig = hilbert(Sig);
end

fft_Sig = fft(Sig)/NN;
fft_Sig = fft_Sig(1: floor(NN/2), :);
ffty_Sig = abs(fft_Sig);

if rowVec == 1
    fft_Sig = fft_Sig.';
    ffty_Sig = ffty_Sig.';
end
