% Example 1 for ACMD

clear, clc, close all
path = '../../common';
addpath(path);
%% Signal
Fs = 1000;
T = 1;
t = 0: 1/Fs: T - 1/Fs;
Sig1 = exp(-0.3*t) .* cos(2*pi * (350*t + 1/2/pi*cos(2*pi*25*t)));
IF1 = 350 - 25*sin(50*pi*t);
Sig2 = exp(-0.6*t) .* cos(2*pi * (250*t + 1/2/pi*cos(2*pi*20*t)));
IF2 = 250 - 20*sin(40*pi*t);

Sig = Sig1 + Sig2;

sigshow(Sig, Fs);
[f, fftSpec] = myfft(Sig, Fs);
fftshow(f, fftSpec);
Spec = mySTFT(Sig, Fs, 512, 32);
stftshow(t, f, Spec);

%% Parameter Setting
tao = 1e-3;
mu = 1e-4;
tol = 1e-8;

%% Component 1 Extraction
[~, findex1] = max(fftSpec);
f1peak = f(findex1);
iniIF1 = f1peak * ones(1, length(Sig));

[IFest1, Sigest1, IAest1] = ACMD(Sig, Fs, iniIF1, tao, mu, tol);

%%


%% Remove Path
rmpath(path);