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

[Sigest1, IFest1, IAest1] = ACMD(Sig, Fs, iniIF1, tao, mu, tol);

%% Component 2 Extraction
SigRes = Sig - Sigest1;
[~, fftSpecRes] = myfft(SigRes, Fs);
[~, findex2] = max(fftSpecRes);
f2peak = f(findex2);
iniIF2 = f2peak * ones(1, length(Sig));

[Sigest2, IFest2, IAest2] = ACMD(SigRes, Fs, iniIF2, tao, mu, tol);

%% Estimated IF
figure(101), hold on, box on;
set(gcf, 'Color', 'w');
plot(t, [IF1; IF2], 'b', 'linewidth', 2);   % true IFs
plot(t, [IFest1; IFest2], 'r--', 'linewidth', 2);   % estimated IF
set(gcf, 'Position', [20 100 415 365]);
xlabel('Time (s)', 'FontSize', 14, 'FontName', 'Times New Roman');
ylabel('Frequency (Hz)', 'FontSize', 14, 'FontName', 'Times New Roman');
axis([0 T 100 Fs/2]);
set(gca, 'Fontsize', 14, 'linewidth', 1, 'fontname', 'Times New Roman');

%% Reconstructed Modes
figure(201)
set(gcf, 'Position', [20 100 500 500], 'Color', 'w');
subplot(211), hold on, box on;
plot(t, Sig1, 'k', 'linewidth', 1);
plot(t, Sigest1, 'b--', 'linewidth', 1);
xlabel('Time (s)', 'FontSize', 14, 'FontName', 'Times New Roman');
ylabel('C1', 'FontSize', 14, 'FontName', 'Times New Roman');
set(gca, 'Fontsize', 14, 'linewidth', 1, 'FontName', 'Times New Roman');

subplot(212), hold on, box on;
plot(t, Sig2, 'k', 'linewidth', 1);
plot(t, Sigest2, 'b--', 'linewidth', 1);
xlabel('Time (s)', 'FontSize', 14, 'FontName', 'Times New Roman');
ylabel('C2', 'FontSize', 14, 'FontName', 'Times New Roman');
set(gca, 'Fontsize', 14, 'linewidth', 1, 'FontName', 'Times New Roman');

%% Remove Path
rmpath(path);