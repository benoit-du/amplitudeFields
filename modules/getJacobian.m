function J = getJacobian(par,X)
%%% 22-11-20        first revision
%%% Benoit Duchet, University of Oxford

%%% calculates the model Jacobian of a two-dimensional WC model

%%% INPUTS
% par:      array of parameter values
% X:        2D state where the Jacobian should be computed

%%% OUTPUTS
% J:        Jacobian at X (2x2 array)

% unpacking model parameters
wIE = par(1);
wEI = par(2);
wEE = par(3);
beta = par(4);
Tau = par(5);
thetaE = par(6);
thetaI = par(7);

% computing the Jacobian at X
f = @(x) 1 / (1 + exp(-beta*(x-1)));
fprime = @(x) beta*f(x)*(1-f(x));
J = 1/Tau*[ -1 + wEE*fprime(thetaE+wEE*X(1)-wIE*X(2)), -wIE*fprime(thetaE+wEE*X(1)-wIE*X(2));...
    wEI*fprime(thetaI+wEI*X(1)), -1];

end