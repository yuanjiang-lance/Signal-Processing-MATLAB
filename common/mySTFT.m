function [Spec, f, t] = mySTFT(Sig, SampFreq, Nfbin, WinLen)
%
% Calculating the Short Time Fourier Transform
%
% ---------- Input ------------
%  Sig: the signal to be analyzed, one row/column
%  SampFreq: sampling frequency (Hz)
%  Nfbin: the number of frequency bins 
%  WinLen: the window length to locate signal in time domain
%
% ---------- Output -----------
%  Spec: the STFT spectrum
%  f: frequency bins
%  t: time bins
%
% Author: Yuan JIANG
% Time: 2024-03-09

if size(Sig, 1) < size(Sig, 2), Sig = Sig.'; end
if (isreal(Sig)), Sig = hilbert(Sig); end

SigLen = length(Sig);
N = 2 * Nfbin - 1;

WinLen = ceil(WinLen / 2) * 2; % ceiling the winidow length into a even number
t = linspace(-1, 1, WinLen)'; % The time of window stands in -1~1 s.
sigma = 0.28;
WinFun = (pi*sigma^2)^(-1/4) * exp((-t.^2)/2/(sigma^2));

Lh = (WinLen - 1)/2; % half of the window length

Spec = zeros(N,SigLen) ;  

for iLoop = 1:SigLen
    
    tau = -min([round(N/2)-1,Lh,iLoop-1]):min([round(N/2)-1,Lh,SigLen-iLoop]);
    temp = floor(iLoop + tau);
    temp1 = floor(Lh+1+tau);
    rSig = Sig(temp);

    rSig = rSig .* conj(WinFun(temp1));
    Spec(1:length(rSig),iLoop) = rSig;
end

Spec = fftshift(fft(Spec),1) / length(rSig); % Calculating the Fourier Transform of each windowed signal and shifting the results.

f = linspace(-SampFreq/2,SampFreq/2, N);
t = (0: SigLen-1)/SampFreq;

Spec = Spec(end-Nfbin+1:end,:);
f = f(end-Nfbin+1:end);
