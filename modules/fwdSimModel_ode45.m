function [X_1,X_2] = fwdSimModel_ode45(par,n,dt,X0,options)
%%% 22-11-20        first revision
%%% Benoit Duchet, University of Oxford

%%% model integration based on ode45, uses the model vector field obtained
%%% from getModelVectorField.m

%%% INPUTS
% par:      array of parameter values
% n:        number of time steps to integrate for
% dt:       integrating time step
% X0:       initial position to start model integration
% options:  ode45 options (see ode45 documentation for more details)

%%% OUTPUTS
% X_1:      trajectory obtained from model integration (first dimension), row vector
% X_2:      trajectory obtained from model integration (second dimension), row vector

t = 0:dt:n*dt;
t = t(1:n);

derivatives = @(t,X) getModelVectorField(X,par);

[~,X] = ode45(derivatives,t,X0,options);
X_1 = X(:,1)';
X_2 = X(:,2)';

end

