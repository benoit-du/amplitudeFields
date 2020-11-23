%%% 22-11-20        first revision
%%% Benoit Duchet, University of Oxford

%%% This script can be used to calculate the isostable amplitude field
%%% associated with the basin of attraction of a two-dimensional dynamical
%%% system with complex conjugate eigenvalues. A two-dimensional
%%% Wilson-Cowan (WC) model is implemented in the relevent functions in the 
%%% "modules" folder as an example. To use this script with a different 
%%% model, you will need to do the following:
%%% - modify the function "getJacobian.m" based on your model
%%% - if you want model integration to be peformed with ode45 (slowler but 
%%%   more accurate), modify the function "getModelVectorField.m"
%%% - if you want model integration to be performed with a custom integration
%%%   scheme, modify "fwdSimModel.m". In the WC example, fwdSimModel calls
%%%   a C implementation of the Euler-Maruyama method (the stochastic part
%%%   is not used for isostable fields). 
%%% - provide the model parameters to used in a matfile (see parFname)
%%% - adapt the script parameters to your problem (in particular the
%%%   coordinates of the grid where the isostable field should be computed, 
%%%   a rough guess of the fixed point position, and the isostable
%%%   computation parameters n  and dt).
%%%
%%% The computation method used in this script is based on 
%%% "Mauroy, A., Mezi?, I., & Moehlis, J. (2013). Isostables, 
%%% isochrons, and Koopman spectrum for the action–angle representation of
%%% stable fixed point dynamics. Physica D: Nonlinear Phenomena, 261,
%%% 19-30.".
%%% The isostable computation parameter n should be chosen large enough
%%% so that the values of the isostable field have converged. Depending on
%%% the numerical accuracy of model forward simulation, large values of n
%%% may lead to numerical instability. 
%%%
%%% Notations
%%% - the state vector is denoted X = [X_1, X_2]
%%% - the isostable amplitude field is denoted by r

clearvars
close all

%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parFname        = 'par.mat'; %name of the model parameter file to load
X_1_grid        = 0.25:1E-3:0.55; %vector of X_1 values where the isostable field should be evaluated 
X_2_grid        = 0.25:1E-3:0.5; %vector of X_2 values where the isostable field should be evaluated 
X_FP_guess      = [0 0]; %rough guess of the fixed point coordinates
dt              = 1E-3; %model integration time step
n               = 5; %number of periods used to approximate the isostable amplitude field (see equation (4) for more details)
useOde45        = false; %if true, model integration will be carried out using fwdSimModel_ode45, if false fwdSimModel will be used (Euler-Maruyama method in the Wilson-Cowan example provided) 
abstol          = 1E-12; %absolute tolerance used by ode45
reltol          = 1E-12; %relative tolerance used by ode45
saveIsoSfield   = true; %a matfile with the obtained isostable amplitude field will be saved if true
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
addpath(['.' filesep 'modules'])

% loading parameters
temp  = load(parFname);
par = temp.par;

%% getting the fixed point and eigenvalues

% locating the fixed point
n_FP = 1E3/dt;
if useOde45
    options = odeset('abstol',abstol,'reltol',reltol);
    [X_1_FP,X_2_FP] = fwdSimModel_ode45(par,n_FP,dt,[X_FP_guess(1),X_FP_guess(2)],options);
else
    [X_1_FP,X_2_FP] = fwdSimModel(par,n_FP,dt,[X_FP_guess(1),X_FP_guess(2)],0);
end
X_star = [X_1_FP(end);X_2_FP(end)];

% obtaining the Jacobian at the fixed point and its eigenvalues
J = getJacobian(par,X_star);
[V,D,W] = eig(J);
eigvals = diag(D);

% checking the eigenvalues are complex conjugates
assert(real(eigvals(1)) ~= 0,'eigenvalues are not complex conjugates: real(eigvals(1)) = 0')
assert(real(eigvals(1)) == real(eigvals(2)),'eigenvalues are not complex conjugates: real(eigvals(1)) ~= real(eigvals(2))')
assert(prod(imag(eigvals)) < 0,'eigenvalues are not complex conjugates: prod(imag(eigvals)) >= 0')
assert(abs(imag(eigvals(1))) == abs(imag(eigvals(2))),'eigenvalues are not complex conjugates: abs(imag(eigvals(1))) ~= abs(imag(eigvals(2)))')

% obtaining the angular frequency and decay
if imag(eigvals(1)) > 0
    i_plus = 1;
else
    i_plus = 2;
end
omega = imag(eigvals(i_plus));
T = 2*pi/omega;
sigma = real(eigvals(i_plus));

%%% used for observables
a = real(V(:,i_plus));
b = -imag(V(:,i_plus));
df1 = [b(2),-b(1)];
df2 = [a(2),-a(1)];

%% computing isostable amplitude field

n_X_1 = length(X_1_grid);
n_X_2 = length(X_2_grid);
r = NaN(n_X_2,n_X_1);
n_dt = round(n*T/dt);

for i_X_1 = 1:n_X_1 
    for i_X_2 = 1:n_X_2
        
        clc
        disp(['i_X_1 = ' num2str(i_X_1) ' / ' num2str(n_X_1)])
        disp(['i_X_2 = ' num2str(i_X_2) ' / ' num2str(n_X_2)])
        
        if useOde45
            [X_1,X_2] = fwdSimModel_ode45(par,n_dt,dt,[X_1_grid(i_X_1),X_2_grid(i_X_2)],options);
        else           
            [X_1,X_2] = fwdSimModel(par,n_dt,dt,[X_1_grid(i_X_1),X_2_grid(i_X_2)],0);
        end
        X = [X_1(n_dt);X_2(n_dt)];
        r(i_X_2,i_X_1) = exp( -sigma*n*T ) * sqrt( dot(df1,X-X_star)^2 + dot(df2,X-X_star)^2 );
    end
end

r = r/abs(dot(df1,a));%normalisation 
tElapsed = toc

%% plotting 

X_1_corners = [X_1_grid(1),X_1_grid(n_X_1)]; %% corners correspond to r(1,1) and r(n_x_2,n_X_1)
X_2_corners = [X_2_grid(1),X_2_grid(n_X_2)];

figure
colormap('jet')
imagesc(X_1_corners,X_2_corners,r)
set(gca,'Ydir','normal')
c = colorbar;
ylabel(c,'$r$ (isostable coordinate)','interpreter','latex');
c.Label.FontSize = 13;
xlabel('$X_1$','interpreter','latex');
ylabel('$X_2$','interpreter','latex');
ftSize = 13;
set(gca,'FontSize',ftSize)

%% saving an isostable amplitude field structure

if saveIsoSfield
    isoS.r = r;
    isoS.X_1_grid = X_1_grid;
    isoS.X_2_grid = X_2_grid;
    isoS.X_1_corners = X_1_corners;
    isoS.X_2_corners = X_2_corners;
    isoS.par = par;
    isoS.usedOde45 = useOde45;
    isoS.n = n;
    isoS.dt = dt;
    isoS.X_star = X_star;
    isoS.sigma = sigma;
    isoS.T = T;
    isoS.tElapsed = tElapsed;
    
    formatOut = 'dd-mmm-yy_HH-MM-ss';
    fname = ['isostableAmpField_' datestr(clock,formatOut)];
    
    matfile saving
    m = matfile(fname,'Writable',true);
    m.isoS = isoS;
end
