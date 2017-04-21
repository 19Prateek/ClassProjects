function [reservoir,numerical,well] = Initialise()

%^^^^^^^^^^^^^^^^^^^^^^^^^^^^USER INPUTS^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
numerical.dt = 0.1;                                  % Time-step size, day
numerical.t0 = 0;                                    % Initial time , days
numerical.t_fin = 500;                               % final time , days
reservoir.Nx = 3;                                    % Number of grids in x
reservoir.Ny = 3;                                    % Number of grids in y
reservoir.muw = 1;                                   % viscosity water, cP
reservoir.muo = 5*ones(reservoir.Nx*reservoir.Ny,1); % viscosity water, cP                             
reservoir.Bw = 1;                                    % FVF water
reservoir.Bo = 1*ones(reservoir.Nx*reservoir.Ny,1);  % FVF oil
reservoir.phi = 0.26;                                % pororsity
reservoir.kx = 1800;                                 % x permeability,mD
reservoir.cw = 2e-6;                                 % water compressibility, psi^-1
reservoir.co = 5e-6;                                 % oil compressibility, psi^-1
reservoir.cr = 3e-6;                                 % rock compressibility, psi^-1
reservoir.Sw0 = 0.2;                                 % Initial Saturation
reservoir.Lx = 1200;                                 % Length in x-dir of reservoir, ft
reservoir.Ly = 600;                                  % Length in y-dir of reservoir, ft
reservoir.h = 200;                                   % Reservoir Thickness, ft
reservoir.z = 0;                                     % Reservoir Depth, ft


%^^^^^^^^^^^^^^^^^^^^^^^^^^WELL INPUTS^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
well.Location = [200 100; 
                 600 300; 
                 1000 500]; 

well.Orientation = [ 0 0 0 ]';              % 1=Horizontal Well,  0=Vertical Well
well.S = [ 0 0 0 ]';                        % skin factor
well.rw = [0.5 0.5 0.5]';                   % well radius, ft

well.nSC = 0;                                        % Number of well schedule changes
well.tSC = 0;                                      % Times of well schedule changes

% Type 1 = constant rate, Type 0 = constant BHP
well.Type = [1 ;
             1 ;
             0 ]; 
         
% Well Rates, ft3/d
well.Rate = [-2000 ;
             3000 ;
             0]*5.61;   
         
% BHP, psi
well.BHP = [0 ;
            0 ;
            800];              

end

