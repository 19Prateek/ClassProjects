function [reservoir,numerical,well] = UserInput

%^^^^^^^^^^^^^^^^^^^^^^^^^^^^USER INPUTS FROM FILE^^^^^^^^^^^^^^^^^^^^^^^^^
reservoir.hM = dlmread('PJ2-Thickness.txt');      % Thickness(from file),ft
reservoir.zM = dlmread('PJ2-Depth.txt');          % Depth(from file), ft
reservoir.phiM = dlmread('PJ2-Porosity.txt');     % Porosity(from file)
reservoir.kM = dlmread('PJ2-Permeability.txt');   % PermX(from file), mD

%^^^^^^^^^^^^^^^^^^^^^^^^^^^^TIME STEP CONTROL^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
numerical.dt = 1;                                 % Time-step(days)
numerical.t0 = 0;                                 % Initial time(days)
numerical.t_fin = 3987;                           % Final time(days)

%^^^^^^^^^^^^^^^^^^^^^^^^^^^^RESERVOIR PROPERTY INPUT^^^^^^^^^^^^^^^^^^^^^^
reservoir.Nx = 80;                                % Number of grids in x
reservoir.Ny = 75;                                % Number of grids in y
reservoir.Lx = 6000;                              % Length in x-dir, ft
reservoir.Ly = 7500;                              % Length in y-dir, ft
reservoir.muw = 0.383211;                         % Water viscosity, cP                             
reservoir.Bw = 1;                                 % FV factor, water
reservoir.rhow = 0.433;                           % Water density, psi/ft
reservoir.rhoo = 53/144;                          % Oil density, psi/ft
reservoir.cw = 3.215e-6;                          % Water comp., psi^-1
reservoir.cr = 1e-6;                              % Rock comp., psi^-1
reservoir.Powoc = 4500;                           % Water pressure at WOC, psi

%^^^^^^^^^^^^^^^^^^^^^^^^^^WELL INPUTS^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% Well locations (x,y), ft
well.Location = [3675 5600;        % Vertical Well 1   
                 3825 3600;        % Vertical Well 1
                 2475 4400;        % Horizontal Well 3
                 2550 4400;        % Horizontal Wekl 3
                 2625 4400;        % Horizontal Well 3
                 2175 2700;        % Horizontal Well 4
                 2250 2700;        % Horizontal Well 4
                 2325 2700;        % Horizontal Well 4
                 1125 1100;        % Vertical Well 5
                 450  3100];       % Vertical Well 6                

% Skin factor(dimensionless constant, -ve for stimulaed wells, +ve for
% damaged wells)
well.S = [0 0 0 0 0 0 0 0 0 0]';

% Well radius, ft
well.rw = [0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25]';  

% Well Orientaion(1 = Horizontal Well, 0 = Vertical Well)
well.Orientation = [0 0 1 1 1 1 1 1 0 0]';        



%^^^^^^^^^^^^^^^^^^^^^^^^^^SCHEDULE DEPENDENT WELL PARAMETERS^^^^^^^^^^^^^^
                    % Each column corresponds to a different schedule

% Well Schedule Change parameters
well.nSC = 1;                            % Number of schedule changes
well.tSC = 1826;                         % Time in days

% Well Type 1 = constant rate and Well Type 0 = constant BHP
well.Type = [1 0;
             1 0;
             1 0;
             1 0;
             1 0;
             1 0;
             1 0;
             1 0;
             0 1;
             0 1]; 
         
% Well BHP(psi) for constant bottom hole pressure wells
well.BHP = [0 100;
            0 75;
            0 50;
            0 50;
            0 50;
            0 100;
            0 100;
            0 100;
            100 0;
            150 0]; 
        
% Well rates(bbl/d) for constant rate wells negative for production, positive for injection        
well.Rate = [-300 0;
             -150 0;
             -400/3 0;
             -400/3 0;
             -400/3 0;
             -250/3 0;
             -250/3 0;
             -250/3 0;
               0 1500;
               0 1500]*5.615;   
           
end

