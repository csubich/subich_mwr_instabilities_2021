function [Dx_phi, Dx_u] = fd_ops(Nx,Lx)

% Make finite-difference operators for the staggered grid
dx = Lx/Nx;


% Helper matrices for defining differential operators
Ptz0 = sparse(1:Nx,1:Nx,ones(Nx,1)); % Identity, evaluate the point 'here'
Ptp1 = circshift(Ptz0,1,2); % Evaluate the point one to the right of 'here'
% Ptp2 = circshift(Ptz0,2,2); % Two to the right
Ptm1 = circshift(Ptz0,-1,2); % One to the left
% Ptm2 = circshift(Ptz0,-2,2); % One to the left
% ZMAT = sparse([],[],[],Nx,Nx); % Zero matrix

% Differential operators
Dx_phi = (Ptp1 - Ptz0)/dx; % Central difference of phi points, which lives at u points
% Avg_phi_u = (Ptp1 + Ptz0)/2; % Average from phi to u grid
Dx_u = (Ptz0 - Ptm1)/dx; % Central difference u points, living at phi points to the left
% Avg_u_phi = (Ptz0 + Ptm1)/2; % Average from u to phi grid