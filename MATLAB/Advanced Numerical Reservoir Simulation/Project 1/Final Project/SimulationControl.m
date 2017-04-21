function [dt,t0,t_fin,nt,rho,zM,kM,phiM,hM,Ny,Nx,mu,Bw] = SimulationControl()

dt = 1;                                  % Time-step size, day
t0 = 0;                                  % Initial time , days
t_fin = 1000;                            % final time , days
nt = round((t_fin-t0)/dt);               % Number of time steps
rho = 0.433;                             % Pressure gradient(water), psi/ft
zM = dlmread('PJ1-Depth.txt');           % Depth Matrix
kM = dlmread('PJ1-Permeability.txt');    % Perm Matrix
phiM = dlmread('PJ1-Porosity.txt');      % Porosity Matrix
hM = dlmread('PJ1-Thickness.txt');       % Thickness Matrix
[Ny,Nx] = size(zM);                      % Number of blocks in x and y dir
mu = 1.95;                               % viscosity, cP              
Bw = 1;                                  % Formation volume factor

end

