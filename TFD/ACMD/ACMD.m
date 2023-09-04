function [Sigest, IFest, IAest] = ACMD(Sig, Fs, iniIF, tao, mu, tol, maxit)
%
% Adaptive Chirp Mode Decomposition (ACMD) without bandwidth adaptation 
%
% Joint-estimation scheme is also integrated in this function, so it could
% extract multiple components simultaneously. If you want to use recursive
% ACMD, just input one initial IF in one execution.
%
% If you want to use ACMD with bandwidth adaptation, please turn to ACMD_adapt.m
%
% ------------- Input ---------------
%  Sig: measured signal, one row/colum vector
%  Fs: sampling frequency (Hz)
%  iniIF: initial instaneous frequencies (IFs), each IF lies in one row
%         ATTENTION: The length of iniIF and Sig must be equal
%  tao: bandwidth controlling parameter, smaller tao results in narrower bandwidth     
%  mu: IF smooth degree controlling parameter, smaller mu results in smoother IF result
%  tol: iteration stopping criterion
%  maxit: maximum iteration number to avoid dead loop, default 300
%
% ------------- Output ---------------
%  Sigest: estimated signal modes, each mode lie in one row
%  IFest: estimated IFs, each IF lies in one row
%  IAest: estimated instantanous amplitudes (IAs), equivalent to the envelope, each IA lies in one row
%
% Notations of each parameter and variable align with the notations in
% original paper, which helps users to better understand the core
% concept of the algorithm.
%
% Author: Yuan JIANG
% Time: 2023-09-02

%% Initialization
if nargin < 7, maxit = 300; end
if length(Sig) ~= size(iniIF, 2)
    error('The length of measured signal and initial IF must be equal!');
end
if size(Sig, 2) > size(Sig, 1), Sig = Sig.'; end

[M, N] = size(iniIF);    % N is the signal length (must equal to the length of IF), M is the mode number
t = (0: N-1) / Fs;
e = ones(N, 1);
e2 = -2 * e;

D = spdiags([e e2 e], 0:2, N-2, N); % 2nd-order difference operator, matrix D in original paper
Ddoub = D'*D;   % matrix D'*D
spzeros = spdiags(zeros(N, 1), 0, N-2, N);
PHIm = [D spzeros; spzeros D];   % matrix PHI in original paper
PHI = kron(eye(M), PHIm);   % Due to the joint-estimation scheme, all PHIm should be integrated in diagonal
PHIdoub = PHI' * PHI;   % matrix PHI'*PHI

IFitset = zeros(M, N, maxit);    % Estimated IFs in each iteration
Sigitset = zeros(M, N, maxit);   % Estimated signal components in each iteration
IAitset = zeros(M, N, maxit);    % Estimated IAs in each iteration

%% Iteration
it = 1;
sDif = tol + 1; % sDif is the energy difference between two consecutive iterations
IF = iniIF;

while (sDif > tol && it <= maxit)
    
    K = sparse(N, 2*N*M);
    for i = 1: M
        cosm = cos(2 * pi * cumtrapz(t, IF(i, :)));
        sinm = sin(2 * pi * cumtrapz(t, IF(i, :)));
        Cm = spdiags(cosm(:), 0, N, N);
        Sm = spdiags(sinm(:), 0, N, N);
        Km = [Cm Sm];
        % Due to the joint-estimation scheme, all Km should be integrated
        % in a row as K = [K1, K2, ... KM]
        startCol = (i-1) * 2*N + 1;
        endCol = startCol + 2*N - 1;
        K(:, startCol: endCol) = Km;
    end
    Kdoub = K' * K;
    
    y = (1/tao * PHIdoub + Kdoub) \ (K' * Sig); 
    
    for i = 1: M
        % updating demodulated sigals
        ym = y((i-1)*2*N+1: i*2*N);
        Km = K(:, (i-1)*2*N+1: i*2*N);
        Sigitset(i, :, it) = Km * ym;
        
        % IF refinement
        alpham = ym(1: N);
        betam = ym(N+1: end);   % two demodulated quadrature signals
        dalpham = Differ(alpham, 1/Fs);
        dbetam  = Differ(betam, 1/Fs);  % derivative of demodulated signals
        dIFm = (betam.*dalpham - alpham.*dbetam) ./ (2*pi* (alpham.^2 + betam.^2));
        dIFm = (1/mu*Ddoub + speye(N)) \ dIFm;
        IF(i, :) = IF(i, :) + dIFm.';
        
        IAitset(i, :, it) = sqrt(alpham.^2 + betam.^2);
    end
    IFitset(:, :, it) = IF;
    
    % convergence criterion
    if it > 1
        sDif = 0;
        for i = 1: M
            sDif = sDif + (norm(Sigitset(i,:,it) - Sigitset(i,:,it-1)) ...
                / norm(Sigitset(i,:,it-1))) .^2;
        end
    end
    it = it + 1;
    
end

it = it - 1;    % final iteration
IFest = IFitset(:, :, it);   % estimated IF
Sigest = Sigitset(:, :, it); % estimated signal component
IAest = IAitset(:, :, it); % estimated IA
