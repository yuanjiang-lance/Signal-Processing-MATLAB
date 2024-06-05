function [Spec, f, t] = mySFFT(sigfft, SampFreq, Nt, WinLen)
% 
% Calculating the Short Frequency Fourier Transform 
% ------- Inpuut --------
%  sigfft : frequency domain signal to be analyzed, one row/column
%  SampFreq : sample Frequency of time domain signal
%  Nt : signal length in time domain
%  WinLen : window length
%
% ------- Output --------
%  Spec : time-frequency spectrum
%  f : frequency bins
%  t : time bins
%
% Author: Yuan JIANG
% Time: 2024-05-26

if size(sigfft, 1) < size(sigfft, 2), sigfft = sigfft.'; end

SigLen = Nt;
N = length(sigfft);           
df = [0:(N-1)]/(N-1)*SampFreq/2;


%%
WinLen = ceil(WinLen / 2) * 2; % ceiling the winidow length into a even number
f = linspace(-1, 1, WinLen)'; % The frequency of window stands in -1~1 Hz.

WinFun = exp(log(0.005) * f.^2 );
WinFun = WinFun / norm(WinFun);
Lh = (WinLen - 1)/2; % half of the window length

fft_len = length(df);
Spec = zeros(SigLen, fft_len);  
conjSpec = Spec;

for iLoop = 1:fft_len
   
    tau = -min([round(fft_len/2)-1,Lh,iLoop-1]):min([round(fft_len/2)-1,Lh,fft_len-iLoop]);  % signal span
    temp = floor(iLoop + tau);
    sSig = sigfft(temp);

    temp1 = floor(Lh+1+tau);    % window span
    sSig = sSig .* conj(WinFun(temp1)); % X(f)* complex conjugate of window
    Spec(1:length(sSig),iLoop) = sSig;  % windowed analytic signal
    conjSpec(:,iLoop)  =  fliplr(Spec(:,iLoop));
end

iSpec  = [conjSpec(1:end-1,:);Spec];
iLen = length(iSpec);
Spec = iSpec(round(iLen/2):end,:);

Spec = ifft(Spec);
Spec = abs(Spec)/2/pi;

[SigLen,nLevel] = size(Spec);

f = (0:nLevel-1) / (nLevel-1) * SampFreq/2;  % frequency in TF plane
t = (0: SigLen-1) / SampFreq;      % time in TF plane
Spec = Spec' / length(sSig);


