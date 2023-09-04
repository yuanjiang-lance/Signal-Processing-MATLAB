function [Sigest, IFest, IAest] = ICCD(Sig, Fs, iniIF, orderIA, lambda, orderIF)
%
% Intrinsic Chirp Component Decomposition (ICCD)
%
% --------------- Input ----------------
%  Sig: measured signal, one row/column vector, real/complex signal
%  Fs: sampling frequency (Hz)
%  iniIF: initial instantaneous frequencies (IFs), each IF lies in one row
%         ATTENTION: The length of iniIF and Sig must be equal
%  orderIA: the order of Fourier series covering instantaneous amplitude (IA). 
%           Higer orderIA results in larger bandwidth
%  lambda: Tikhonov regularization parameter. 
%          Higer lambda results in insensitivity of noise
%  orderIF (optional): the order of Fourier series for IF fitting. 
%                      If orderIF is omitted, iniIF should be smoothed IFs with
%                      high accuracy (e.g. after fitting or low-pass
%                      filtering). In this case, the output IFest will
%                      equal to iniIF.
%                      
% --------------- Output -----------------
%  Sigest: estimated signal modes, each mode lie in one row
%  IFest: estimated IFs, each IF lies in one row
%  IAest: estimated instantanous amplitudes (IAs), equivalent to the envelope, each IA lies in one row
%
% Author: Yuan JIANG
% Time: 2023-09-03

%% Initialization
if nargin < 6, fineIF = true; else, fineIF = false; end
if length(Sig) ~= size(iniIF, 2)
    error('The length of measured signal and initial IF must be equal!');
end
if size(Sig, 2) > size(Sig, 1), Sig = Sig.'; end
if isreal(Sig), realSig = true; else, realSig = false; end

[M, N] = size(iniIF);    % N is the signal length (must equal to the length of IF), M is the mode number
t = (0: N-1) / Fs;

if fineIF
    IFest = iniIF;
    phase = 2*pi*cumtrapz(t, IFest, 2);
else
    [IFest, phase] = IFfit_overFourier(iniIF, Fs, orderIF);
end

%% Constructing Fourier Matrix
f0 = Fs / (2*N);    % base frequency
K = zeros(N, 2*orderIA + 1);    % K matrix in original paper
for i = 1: orderIA + 1
    K(:, i) = cos(2*pi * (i-1) * f0 * t);
end
for i = orderIA + 2: 2*orderIA + 1
    K(:, i) = sin(2*pi * (i-orderIA-1) * f0 * t);
end

if realSig
    H = zeros(N, (2*orderIA+1)*2*M);
    for i = 1: M
        Ci = cos(phase(i, :));
        H(:, (i-1)*2*(2*orderIA+1)+1: (2*i-1)*(2*orderIA+1)) = spdiags(Ci(:), 0, N, N) * K;
        Si = sin(phase(i, :));
        H(:, (2*i-1)*(2*orderIA+1)+1: 2*i*(2*orderIA+1)) = spdiags(Si(:), 0, N, N) * K;
    end
else
    H = zeros(N, (2*orderIA+1)*M);
    for i = 1: M
        Ci = exp(1j*phase(i, :));
        H(:, (i-1)*(2*orderIA+1)+1: i*(2*orderIA+1)) = spdiags(Ci(:), 0, N, N) * K;
    end
end

%% Recovering signal modes
I = speye(size(H, 2));
y = (H'*H + lambda*I) \ (H'*Sig);

Sigest = zeros(M, N);
IAest = zeros(M, N);

for i = 1: M
    if realSig
        Sigest(i, :) = H(:, (i-1)*2*(2*orderIA+1)+1: 2*i*(2*orderIA+1)) * ...
            y((i-1)*2*(2*orderIA+1)+1: 2*i*(2*orderIA+1));
        a = K * y((i-1)*2*(2*orderIA+1)+1: (2*i-1)*(2*orderIA+1));
        b = K * y((2*i-1)*(2*orderIA+1)+1: 2*i*(2*orderIA+1));
        IAest(i, :) = sqrt(a.^2 + b.^2);
    else
        Sigest(i, :) = H(:, (i-1)*(2*orderIA+1)+1: i*(2*orderIA+1)) * ...
            y((i-1)*(2*orderIA+1)+1: i*(2*orderIA+1));
        IAest(i, :) = K * y((i-1)*(2*orderIA+1)+1: i*(2*orderIA+1));
    end
end

end

%%
function [IFfit, phase] = IFfit_overFourier(iniIF, Fs, orderIF, lambda)

% Fitting instantaneous frequencies (IFs) using over Fourier series

% ------------------- Input -----------------------
%  iniIF: initial instantaneous frequencies (IFs), each IF lies in one row
%  Fs: sampling frequency (Hz)
%  orderIF: the order of Fourier series, one number or a vector with the
%           same length of iniIF
%  lambda: Tikhonov regularization parameter, default 5e-5
%
% ------------------- Output ----------------------
%  IFfit: IF after fitting
%  phase: phase corresponding to fitted IF
%
% Author: Yuan JIANG
% Time: 2023-09-03

if nargin < 4, lambda = 5e-5; end
[M, N] = size(iniIF);
if length(orderIF) ~= N && length(orderIF) ~= 1
    error('The length of orderIF must be 1 or equal to the number of IFs.');
end
if length(orderIF) == 1
    orderIF = orderIF * ones(1, M);
end

t = (0: N-1) / Fs;
f0 = Fs / (2*N);    % base frequency


end