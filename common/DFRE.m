function index = DFRE(Spec, delta)
%
% Dual Fast Ridge Estimation algorithm for time-frequency ridge detection
% Only one ridge curve could be extracted under one excecution
% Suitable for time-related signals, e.g. harmonics
% For frequency-related (dispersive) signals, it is recommended to use DFRE_F
%
% ------------ Input -------------
%  Spec: Time-frequency Spectrum
%  delta: searching scope
%
% ------------ Output ------------
%  index: index of frequency bins corresponding to ridge curve
%
% Author: Yuan JIANG
% Time: 2023-08-31

Spec = abs(Spec);
[M, N] = size(Spec);
index = zeros(1, N);
[fmax, tmax] = find(Spec == max(Spec(:)));
fmax = fmax(1);
tmax = tmax(1);
index(tmax) = fmax;

% Extracting right-side ridge points
f0 = fmax;
for j = min(tmax+1, N): N
    low = max(1, f0-delta);
    up  = min(M, f0+delta);
    [~, f0] = max(Spec(low:up, j));
    f0 = f0 + low - 1;
    index(j) = f0;
end

% Extracting left-size ridge points
f1 = fmax;
for j = max(1, tmax-1): -1: 1
    low = max(1, f1-delta);
    up  = min(M, f1+delta);
    [~, f1] = max(Spec(low:up, j));
    f1 = f1 + low - 1;
    index(j) = f1;
end