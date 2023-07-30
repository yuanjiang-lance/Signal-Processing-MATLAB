function ybar = Differ(y,delta)
%
% Computing the derivative of a discrete time series y
%
% ---------- Input -----------
%  y: time series (e.g. a signal)
%  delta: sampling time interval of y
%
% ---------- Output -----------
%  ybar: derivative of y
%
% Author: Yuan JIANG
% Time: 2023-07-30

L = length(y);
ybar = zeros(1, L-2);
for i = 2: L-1
  ybar(i-1) = (y(i+1)-y(i-1)) / (2*delta);
end
ybar = [(y(2)-y(1))/delta, ybar, (y(end)-y(end-1))/delta];
