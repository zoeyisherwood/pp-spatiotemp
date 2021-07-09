function probs = one_over_x(params,x)

% Design matrix for f(x) = a*x^b.
% 
%
% % CC 6.3.01

offSetParam = 100;

probs = params(1).*(x.^(params(2)))+offSetParam;
