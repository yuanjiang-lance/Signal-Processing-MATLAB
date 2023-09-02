function [IFest, Sigest, IAest, taorec] = ACMD_adapt(Sig, Fs, iniIF, tao0, mu, tol, maxit)
%
% Adaptive Chipr Mode Decomposition (ACMD) with bandwidth adaptation
%
% ------------- Input ---------------
%  Sig: measured signal, one row/colum vector
%  Fs: sampling frequency (Hz)
%  iniIF: initial instaneous frequency (IF), one row/column vector
%         ATTENTION: The length of iniIF and Sig must be equal
%  tao0: initial bandwidth controlling parameter, smaller tao0 results in narrower bandwidth     
%  mu: IF smooth degree controlling parameter, smaller mu results in smoother IF result
%  tol: iteration stopping criterion
%  maxit: maximum iteration number to avoid dead loop, default 300
%
% ------------- Output ---------------
%  IFest: estimated IF
%  Sigest: estimated signal mode
%  IAest: estimated instantanous amplitude (IA), equivalent to the envelope
%  taorec: the recording of tao (bandwidth controlling parameter) in each iteration
%
% Notations of each parameter and variable align with the notations in
% original paper, which helps users to better understand the core
% concept of the algorithm.
%
% Author: Yuan JIANG
% Time: 2023-09-01

%% Initialization
if narin < 7, maxit = 300; end
if length(Sig) ~= length(iniIF)
    error('The length of measured signal and initial IF must be equal!');
end
if size(Sig, 1) > size(Sig, 2), Sig = Sig.'; end
if size(iniIF, 1) > size(iniIF, 2), iniIF = iniIF.'; end

N = length(Sig);    % signal length (must equal to length(iniIF))
t = (0: N-1) / Fs;
e = ones(N, 1);
e2 = -2 * e;

D = spdiags([e e2 e], 0:2, N-2, N); % 2nd-order difference operator, matrix D in original paper
Ddoub = D'*D;   % matrix D'*D
spzeros = spdiags(zeros(N, 1), 0, N-2, N);
PHI = [D spzeros; spzeros D];   % matrix PHI in original paper
PHIdoub = PHI' * PHI;   % matrix PHI'*PHI

IFitset = zeros(maxit, N);    % Estimated IF in each iteration
Sigitset = zeros(maxit, N);   % Estimated signal component in each iteration
taorec = zeros(1, maxit);

%% Iteration
it = 1;
sDif = tol + 1; % sDif is the energy difference between two consecutive iterations
IF = iniIF;
tao = tao0;

while (sDif > tol && it <= maxit)
    
    cosm = cos(2 * pi * cumtrapz(t, IF));
    sinm = sin(2 * pi * cumtrapz(t, IF));
    Cm = spdiags(cosm(:), 0, N, N);
    Sm = spdiags(sinm(:), 0, N, N);
    Km = [Cm Sm];
    Kmdoub = Km' * Km;
    
    % updating demodulated sigals
    ym = (1/tao * PHIdoub + Kmdoub) \ (Km' * Sig(:)); 
    Sigm = Km * ym;
    Sigitset(it, :) = Sigm;
    
    % IF refinement
    alpham  = ym(1: N);
    betam   = ym(N+1: end);     % two demodulated quadrature signals
    dalpham = Differ(alpham, 1/Fs);
    dbetam  = Differ(betam, 1/Fs);  % derivative of demodulated signals
    dIF = (betam.*dalpham - alpham.*dbetam) ./ (2*pi* (alpham.^2 + betam.^2));
    dIF = (1/mu*Ddoub + speye(N)) \ dIF;
    IF = IF + dIF.';
    IFitset(it, :) = IF;
    
    % bandwith adaptation
    tao = tao * (Sigm'*Sigm) / (Sigm'*Sig(:));
    taorec(it) = tao;
    
    % convergence criterion
    if it > 1
        sDif = (norm(Sigitset(it,:) - Sigitset(it-1,:)) / norm(Sigitset(it-1,:))).^2;
    end
    it = it + 1;
    
end

it = it - 1;    % final iteration
IFest = IFitset(it, :);   % estimated IF
Sigest = Sigitset(it, :); % estimated signal component
IAest = sqrt(alpham.^2 + betam.^2).'; % estimated IA
taorec = taorec(1: it);