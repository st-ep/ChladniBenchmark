%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate 1100 Forcing Samples for Chladni Plate (Faster Version)
%   - Precompute everything that does NOT depend on alpha(n,m).
%   - Plot the first two samples immediately.
%   - Split into 1000 Training Samples + 100 Testing Samples.
%   - Save to ChladniData.mat, now including S(x,y) arrays.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc;

%% 1) Basic Setup
L = 8.75 * 0.0254;   % Dimensions in meters
M = 8.75 * 0.0254;
omega   = 55 * pi / M;  % Frequency
t_fixed = 4;            % Time at which to evaluate the solution

gamma = 0.02;  % damping_adjustment
v     = 0.5;

numPoints = 100;
x = linspace(0, L, numPoints);
y = linspace(0, M, numPoints);

n_range = 20;  
m_range = 20;

N_total = 1100;   % total number of samples
N_train = 1000;   % number of training samples
N_test  = 100;    % number of testing samples

% We'll store:
%   alpha_full : (n_range, m_range, N_total)
%   S_full     : (numPoints, numPoints, N_total)  <-- the forcing
%   Z_full     : (numPoints, numPoints, N_total)  <-- the displacement
alpha_full = zeros(n_range, m_range, N_total);
S_full     = zeros(numPoints, numPoints, N_total);
Z_full     = zeros(numPoints, numPoints, N_total);

%% 2) Precompute Terms That Don't Depend on alpha
% ------------------------------------------------

% (A) Wave numbers mu(n), lamda(m)
mu_vals   = ( (1:n_range)' * pi ) / L;   % column vector of length n_range
lam_vals  = ( (1:m_range)' * pi ) / M;   % column vector of length m_range

% (B) cosX(n,i) = cos( mu_vals(n) * x(i) )
cosX = zeros(n_range, numPoints);
for n = 1:n_range
    for i = 1:numPoints
        cosX(n,i) = cos(mu_vals(n) * x(i));
    end
end

% (C) cosY(m,j) = cos( lam_vals(m) * y(j) )
cosY = zeros(m_range, numPoints);
for m = 1:m_range
    for j = 1:numPoints
        cosY(m,j) = cos(lam_vals(m) * y(j));
    end
end

% (D) centerFactor(n,m) = cos(mu_n*(L/2)) * cos(lam_m*(M/2))
centerFactor = zeros(n_range, m_range);
for n = 1:n_range
    for m = 1:m_range
        centerFactor(n,m) = cos(mu_vals(n)*(L/2)) * cos(lam_vals(m)*(M/2));
    end
end

% (E) beta(n,m) = sqrt(mu^2 + lam^2 + 3*v^2 - gamma^4)
beta_nm = zeros(n_range, m_range);
for n = 1:n_range
    for m = 1:m_range
        beta_nm(n,m) = sqrt( mu_vals(n)^2 + lam_vals(m)^2 + 3*v^2 - gamma^4 );
    end
end

% (F) timeInt(n,m) = \int_0^t sin(omega*(tau - t)) * exp(-gamma^2 + v^2*tau)
%                                  * sin(beta_nm(n,m)*tau) d\tau
timeInt = zeros(n_range, m_range);
for n = 1:n_range
    for m = 1:m_range
        
        current_beta = beta_nm(n,m);
        
        integrand = @(tau) ...
            sin( omega*(tau - t_fixed) ) .* ...
            exp( -gamma^2 + v^2 * tau ) .* ...
            sin( current_beta * tau );
        
        timeInt(n,m) = integral(integrand, 0, t_fixed);
    end
end

% (G) modeFactor(n,m) = (v^2 / beta_nm(n,m)) * timeInt(n,m) * (4/(L*M)) * centerFactor(n,m)
modeFactor = zeros(n_range, m_range);
for n = 1:n_range
    for m = 1:m_range
        modeFactor(n,m) = ...
            (v^2 / beta_nm(n,m)) * timeInt(n,m) * (4/(L*M)) * centerFactor(n,m);
    end
end

%% 3) Main Loop: Generate Data
%    - For each sample, we create alpha_k, then:
%      1) Build S_k(i,j) = sum alpha_k(n,m)*cosX(n,i)*cosY(m,j)
%      2) Build Z_k(i,j) = sum alpha_k(n,m)*cosX(n,i)*cosY(m,j)*modeFactor(n,m)
% ------------------------------------------------------------------------

for k = 1:N_total
    % (A) Generate random alpha-coefficients
    alpha_k = 0.01 * randn(n_range, m_range);
    alpha_full(:,:,k) = alpha_k;
    
    % (B) Compute S(i,j) for each grid point
    %        S_k = sum_{n,m} alpha_k(n,m)*cosX(n,i)*cosY(m,j)
    S_k = zeros(numPoints, numPoints);
    for i = 1:numPoints
        for j = 1:numPoints
            temp_sum = 0;
            for n = 1:n_range
                cx = cosX(n,i);
                for m = 1:m_range
                    temp_sum = temp_sum + alpha_k(n,m)*cx*cosY(m,j);
                end
            end
            S_k(i,j) = temp_sum;
        end
    end
    S_full(:,:,k) = S_k;
    
    % (C) Compute Z(i,j) for each grid point
    %        Z_k = sum_{n,m} alpha_k(n,m)*cosX(n,i)*cosY(m,j)*modeFactor(n,m)
    Z_k = zeros(numPoints, numPoints);
    for i = 1:numPoints
        for j = 1:numPoints
            temp_sum = 0;
            for n = 1:n_range
                cx = cosX(n,i);
                for m = 1:m_range
                    temp_sum = temp_sum + alpha_k(n,m)*cx*cosY(m,j)*modeFactor(n,m);
                end
            end
            Z_k(i,j) = temp_sum;
        end
    end
    Z_full(:,:,k) = Z_k;
    
    % (D) Plot the first few samples as soon as they're computed
    if k <= 4
        figure('Color','black');
        contour(x, y, Z_k, 'LevelList', 0, 'LineColor','white','LineWidth',3);
        xlabel('X axis','Color','white');
        ylabel('Y axis','Color','white');
        title(['Sample #' num2str(k) ', \omega = ' num2str(omega) ' Hz'],'Color','white');
        set(gca,'Color','black','XColor','white','YColor','white');
        
        drawnow;  % show plot immediately
    end
end

%% 4) Split into Training and Testing, then Save
% ----------------------------------------------
alpha_train = alpha_full(:,:,1:N_train);
alpha_test  = alpha_full(:,:,N_train+1:end);

% The forcing S is 3D: (numPoints x numPoints x N_total).
S_train = S_full(:,:,1:N_train);
S_test  = S_full(:,:,N_train+1:end);

Z_train = Z_full(:,:,1:N_train);
Z_test  = Z_full(:,:,N_train+1:end);

save('ChladniData.mat', ...
    'alpha_train','alpha_test', ...
    'S_train','S_test', ...
    'Z_train','Z_test', ...
    'x','y', ...
    '-v7.3');
