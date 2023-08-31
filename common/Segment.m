function SigSegment = Segment(Sig, Fs, TStart, TEnd)
%
% Truncating a segment from a long signal based on time
%
% -------- Input ----------
%  Sig: original signal, one row/column
%  Fs: sampling frequency (Hz)
%  TStart: start time of truncation (s)
%  TEnd: end time of truncation (s)
%
% -------- Output ---------
%  SigSegment: signal segment, one column
%
% Author: Yuan JIANG
% Time: 2023-08-31

if size(Sig, 1) == 1
    Sig = Sig.';
end
t = (0: size(Sig,1)-1) / Fs;
[~, NStart] = min(abs(t-TStart)); % regarding the closest time point to TStart as the starting point
N = ceil((TEnd - TStart)*Fs); % ceiling the total length of signal segment
NEnd = NStart + N - 1;
SigSegment = Sig(NStart:NEnd, :);