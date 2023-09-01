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
% Time: 2023-08-31

L = length(y);
ybar = (y(3:end) - y(1:end-2)) / (2*delta);
ybar = [(y(2)-y(1))/delta, ybar, (y(end)-y(end-1))/delta];
