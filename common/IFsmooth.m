function newIF = IFsmooth(IF, mu)
%
% Curve smooth for instantaneous frequencies (IFs)
%
% ------------- Input ---------------
%  IF: instantaneous frequencies, each IF must be listed in one row
%  mu: smooth degree controlling parameter, smaller beta results in smoother output IF
%
% ------------- Output --------------
%  newIF: smoothed IF
%
% Author: Yuan JIANG
% Time: 2024-05-30

[M, N] = size(IF);  % M is the number of IFs (components), N is the length of each IF (component)
if length(mu) == 1
    mu = mu * ones(1, M);
elseif length(mu) ~= M
    error('The length of smoothness control parameter should equal the row number of IF');
end
e = ones(N, 1);
e2 = -2 * e;
PHI = spdiags([e e2 e], 0:2, N-2, N);   % 2nd-order difference operator
PHIdoub = PHI' * PHI;
newIF = zeros(M, N);
for i = 1: M
    newIF(i, :) = (2/mu(i) * PHIdoub + speye(N)) \ IF(i, :).';
end
