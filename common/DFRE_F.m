function index = DFRE_F(Spec, delta)
%
% Dual Fast Ridge Estimation algorithm for time-frequency ridge detection
% Only one ridge curve could be extracted under one excecution
% Suitable for frequency-related (dispersive) signals, e.g. Lamb waves
% For time-related signals, it is recommended to use DFRE
%
% ------------ Input -------------
%  Spec: Time-frequency Spectrum
%  delta: searching scope
%
% ------------ Output ------------
%  index: index of time bins corresponding to ridge curve
%
% Author: Yuan JIANG
% Time: 2023-09-04

Spec = abs(Spec);
[M, N] = size(Spec);
index = zeros(1, M);
[fmax,tmax] = find(Spec == max(Spec(:)));
fmax = fmax(1);
tmax = tmax(1);
index(fmax) = tmax;

% Extracting ridge points in higher frequency region
t0 = tmax;
for i = (min(fmax+1, M)): M
    low = max(1, t0-delta);
    up = min(N, t0+delta);
    [~, t0] = max(Spec(i, low:up));
    t0 = t0 + low - 1;
    index(i) = t0;
end

t1 = tmax;
for i = (max(1, fmax-1)): -1: 1
    low = max(1, t1-delta);
    up = min(N, t1+delta);
    [~, t1] = max(Spec(i, low:up));
    t1 = t1 + low - 1;
    index(i) = t1;
end