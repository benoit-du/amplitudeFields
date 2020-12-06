%%% 22-11-20        first revision
%%% Benoit Duchet, University of Oxford

%%% This script can be used to calculate the Hilbert amplitude field
%%% associated with the first dimension of a two-dimensional dynamical
%%% system. The present code assumes complex conjugate eigenvalues, 
%%% however this requirement could be lifted with simple modifications 
%%% to the code). A two-dimensional Wilson-Cowan (WC) model is implemented 
%%% in the relevent functions in the "modules" folder as an example. To
%%% use this script with a different model, you will need to do the 
%%% following:
%%% - modify the function "getJacobian.m" based on your model
%%% - modify "fwdSimModel.m". In the WC example, fwdSimModel calls
%%%   a C implementation of the Euler-Maruyama method. This is faster 
%%%   than Matlab code.
%%% - provide in a matfile the model parameters to use (see parFname), 
%%%   including the noise standard deviation as the last element of par 
%%% - adapt the script parameters to your problem (in particular the
%%%   coordinates of the grid where the isostable field should be computed, 
%%%   a rough guess of the fixed point position, as well as nTraj and 
%%%   nPeriods).
%%%
%%% Notations
%%% - the state vector is denoted X = [X_1, X_2]
%%% - the Hilbert amplitude field is called HbAmpField


clearvars
close all

%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parFname            = 'par.mat'; %name of the model parameter file to load
X_1_grid            = 0.25:1E-3:0.55; %vector of X_1 values where the isostable field should be evaluated (must be regularly spaced)
X_2_grid            = 0.25:1E-3:0.5; %vector of X_2 values where the isostable field should be evaluated (must be regularly spaced)
X_FP_guess          = [0 0]; %rough guess of the fixed point coordinates
dt                  = 1E-3; %model integration time step
useMeanCentering    = true; %if true: X_1 is centered by substracting its mean, if false: X_1 is centered by substracting the first coordinate of the fixed point
nTraj               = 30; %number of trajectories generated for averaging
nPeriods            = 35; %number of periods per trajectories
clipEnds            = true; %to clip edge artefacts in Hilbert amplitude time series
pctClip             = 0.5; %how much to clip, in percent of the length of a single trajectory
saveHbField         = true; %a matfile with the obtained Hilbert amplitude field will be saved if true
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
addpath(['.' filesep 'modules'])

% loading parameters
temp  = load(parFname);
par = temp.par;

%% getting the fixed point and eigenvalues

% locating the fixed point
n_FP = 1E3/dt;
[X_1_FP,X_2_FP] = fwdSimModel(par,n_FP,dt,[X_FP_guess(1),X_FP_guess(2)],0);
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

% obtaining the angular frequency
if imag(eigvals(1)) > 0
    i_plus = 1;
else
    i_plus = 2;
end
omega = imag(eigvals(i_plus));
T = 2*pi/omega;


%% generating trajectories

% sampling boundaries
X_1_min = X_1_grid(1);
X_1_max = X_1_grid(end);
X_2_min = X_2_grid(1);
X_2_max = X_2_grid(end);

% initialising arrays
nMax = round(nPeriods*T/dt);%number of time steps to integrate trajectories for
if clipEnds
    nClip = ceil(nMax*pctClip/100);
    nClipped = nMax - 2*nClip;
    X_1_traj = NaN(nTraj,nClipped);
    hAmp = NaN(nTraj,nClipped);
    X_2_traj = NaN(nTraj,nClipped);
else
    X_1_traj = NaN(nTraj,nMax);
    hAmp = NaN(nTraj,nMax);
    X_2_traj = NaN(nTraj,nMax);
end

disp('generating traj **********************')
disp(newline)

for k = 1:nTraj
    progressUpdate(k,nTraj)
    
    % initial position sampling
    X_1_0 = uniRand(X_1_min,X_1_max);
    X_2_0 = uniRand(X_2_min,X_2_max);
    
    % model integration
    [X_1_traj_temp,X_2_traj_temp] = fwdSimModel(par,nMax,dt,[X_1_0,X_2_0],par(end));
    
    % centering the X_1 component and obtaining its hilbert amplitude
    if useMeanCentering
        X_1_centered = X_1_traj_temp-mean(X_1_traj_temp);
    else
        X_1_centered = X_1_traj_temp-X_1_FP;
    end
    hTemp = abs(hilbert(X_1_centered));
    
    % clipping ends to remove edge artefacts
    if clipEnds
        clipIdx = nClip+1:nMax-nClip;
        X_1_traj(k,:) = X_1_traj_temp(clipIdx);
        X_2_traj(k,:) = X_2_traj_temp(clipIdx);
        hAmp(k,:) = hTemp(clipIdx);
    else
        X_1_traj(k,:) = X_1_traj_temp;
        X_2_traj(k,:) = X_2_traj_temp;   
        hAmp(k,:) = hTemp;
    end    
end

disp(newline)


%% bining of trajectory points and averaging Hilbert amplitude in each bin

% creating bin vectors XX_1 and XX_2
X_1_min = X_1_grid(1);
X_1_max = X_1_grid(end);
X_2_min = X_2_grid(1);
X_2_max = X_2_grid(end);
nX_1 = length(X_1_grid);
nX_2 = length(X_2_grid);
dX_1 = (X_1_max - X_1_min) / (nX_1-1);
XX_1 = linspace(X_1_min-dX_1/2, X_1_max+dX_1/2, nX_1+1);
dX_2 = (X_2_max - X_2_min) / (nX_2-1);
XX_2 = linspace(X_2_min-dX_2/2, X_2_max+dX_2/2, nX_2+1);

HbAmpAll = cell(nX_2,nX_1);
smoothFact = 100;
n_smooth = round(T/dt/smoothFact);

disp(newline)
disp('averaging **********************')
disp(newline)
for k = 1:nTraj
    progressUpdate(k,nTraj)
    
    %smoothing
    hAmp_sm = movmean(hAmp(k,:),n_smooth);
    X_1_sm = movmean(X_1_traj(k,:),n_smooth);
    X_2_sm = movmean(X_2_traj(k,:),n_smooth);
    
    %binning 
    for i = 1:nX_2
        idx_i = XX_2(i)<X_2_sm & X_2_sm<XX_2(i+1);
        for j = 1:nX_1
            idx_j = XX_1(j)<X_1_sm & X_1_sm<XX_1(j+1);
            idx_ij = idx_i & idx_j;
            if sum(idx_ij)>1
                HbAmpAll{i,j} = [HbAmpAll{i,j},hAmp_sm(idx_ij)];
            end
        end
    end
end

% averaging in each bin
HbAmpField = cellfun(@mean,HbAmpAll);

% grid vectors (bin center points), should be the same as user defined grid vectors
X_1_centered = (XX_1(1:end-1)+XX_1(2:end))/2;
X_2_centered = (XX_2(1:end-1)+XX_2(2:end))/2;

disp(newline)
tElapsed = toc


%% plotting

X_1_corners = [min(X_1_centered) max(X_1_centered)];
X_2_corners = [min(X_2_centered) max(X_2_centered)];

figure
colormap(jet); 
imagesc(X_1_corners,X_2_corners,HbAmpField,'AlphaData',~isnan(HbAmpField)) 
set(gca,'YDir','normal');
c = colorbar;
ylabel(c,'Hilbert amplitude','interpreter','latex');
c.Label.FontSize = 13;
xlabel('$X_1$','interpreter','latex');
ylabel('$X_2$','interpreter','latex');
ftSize = 13;
set(gca,'FontSize',ftSize)


%% saving a Hilbert amplitude field structure

if saveHbField
    hilbertField.hAmpField = HbAmpField;
    hilbertField.X_1_grid = X_1_centered;
    hilbertField.X_2_grid = X_2_centered;
    hilbertField.X_1_corners = X_1_corners;
    hilbertField.X_2_corners = X_2_corners;
    hilbertField.dt = dt;
    hilbertField.nTraj = nTraj;
    hilbertField.nPeriods = nPeriods;
    hilbertField.useMeanCentering = useMeanCentering;
   
    formatOut = 'dd-mmm-yy_HH-MM-ss';
    fname = ['HilbertAmpField_' datestr(clock,formatOut)];
    save(fname,'hilbertField') 
end