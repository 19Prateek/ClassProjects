function [reservoir,numerical,well] = InitialiseBL()


%^^^^^^^^^^^^^^^^^^^^^^^^^^User reservoir inputs ^^^^^^^^^^^^^^^^^^^^^^^^^^
numerical.dt = 0.1;                                % Time-step size, day
numerical.t0 = 0;                                  % Initial time , days
numerical.t_fin = 190;                             % final time , days
reservoir.Nx = 102;                                % Number of x grids , days
reservoir.Ny = 1;                                  % Number of y grids , days
reservoir.muw = 0.383211;                          % viscosity water, cP
reservoir.muo = 1.03*ones(reservoir.Nx*reservoir.Ny,1);   % viscosity oil, cP                           
reservoir.Bw = 1;                                  % FVF water
reservoir.Bo = 1*ones(reservoir.Nx*reservoir.Ny,1);% FVF oil
reservoir.rhow = 0;                                % density water 
reservoir.rhoo = 0;                                % density oil
reservoir.phi = 0.25;                              % porosity
reservoir.kx = 50;                                 % permeability(x)
reservoir.cw = 0;                                  % water compressibility, psi^-1
reservoir.co = 0;                                  % oil compressibility, psi^-1
reservoir.cr = 1e-6;                               % total compressibility, psi^-1
reservoir.Sw0 = 0.2;                               % Intial water satruation
reservoir.Lx = 510;                                % Length in x-dir of reservoir, ft
reservoir.Ly = 50;                                 % Length in y-dir of reservoir, ft
reservoir.h = 25;                                  % reservoir thickness(ft)
reservoir.z = 0;                                   % reservoir depths (ft)
reservoir.hM = 25*ones(reservoir.Nx,reservoir.Ny);      % Thickness(from file),ft
reservoir.zM = 0*ones(reservoir.Nx,reservoir.Ny);         % Depth(from file), ft
reservoir.phiM =  0.25*ones(reservoir.Nx,reservoir.Ny);     % Porosity(from file)
reservoir.kM =  50*ones(reservoir.Nx,reservoir.Ny);   % PermX(from file), mD

%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^Wells Module^^^^^^^^^^^^^^^^^^^^^^^^^^^
well.Location = [0.5 25; 509.5 25];                % [2.5 25 ; 507.5 25]
well.Orientation = [0 0]';                         % 1=Horizontal Well,  0=Vertical Well
well.S = [0 0]';                                   % skin factor
well.rw = [0.25 0.25]';                            % well radius, ft
well.nSC = 0;                                      % Number of well schedule changes
well.tSC = 0;                                      % Times of well schedule changes
well.Type = [1;1];                                 % Type 1 = constant rate, Type 0 = constant BHP           
well.Rate = [50;-50]*5.615;                        % Well Rates, ft3/d               
well.BHP = [0;0];                                  % BHP, psi
                   
end

