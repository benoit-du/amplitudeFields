function [ Xdot ] = getModelVectorField(X,par)
%%% 22-11-20        first revision
%%% Benoit Duchet, University of Oxford

%%% calculated the vector field of a two-dimensional WC model

%%% INPUTS
% par:      array of parameter values
% X:        2D state where the vector field should be evaluated

%%% OUTPUTS
% Xdot:     2D vector field at X

% unpacking model parameters
wIE = par(1);
wEI = par(2);
wEE = par(3);
beta = par(4);
Tau = par(5);
thetaE = par(6);
thetaI = par(7);

% calculating the vector field at X
E = X(1);
I = X(2);
Edot = (-E+1/(1+exp(-beta*(thetaE-wIE*I+wEE*E-1))))/Tau;
Idot = (-I+1/(1+exp(-beta*(thetaI+wEI*E-1))))/Tau;
Xdot = [Edot;Idot];

end

