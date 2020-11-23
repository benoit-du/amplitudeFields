function [X_1,X_2] = fwdSimModel(par,n,dt,X0,zeta)
%%% 22-11-20        first revision
%%% Benoit Duchet, University of Oxford

%%% wrapper for mex function implementing model forward simulation with
%%% a Euler-Maruyama scheme or equivalent. This function and the mex
%%% function are written for a two-dimensional WC model.

%%% INPUTS
% par:      array of parameter values
% n:        number of time steps to integrate for
% dt:       integrating time step
% X0:       initial position to start model integration
% zeta:  	noise standard deviation

%%% OUTPUTS
% X_1:      trajectory obtained from model integration (first dimension)
% X_2:      trajectory obtained from model integration (second dimension)

% unpacking parameters
wIE = par(1);
wEI = par(2);
wEE = par(3);
beta = par(4);
Tau = par(5);
thetaE = par(6);
thetaI = par(7);

% model forward simulation
[X_1,X_2] = fwdSimEulerMaruyama_WC(wIE,wEI,wEE,beta,Tau,thetaE,thetaI,n,zeta,dt,X0(1),X0(2));

end