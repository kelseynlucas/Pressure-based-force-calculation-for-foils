function args = parameters2()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters2.m  Function to provide parameters to queen2
%
% Parameter description:
% datafolder: directory containing all data (input/output) files
%
% inroot: input velocity filename root string (i.e. filename string
%         that precedes file number)
%
% outroot: output pressure filename 
%
% blanking: set to 1 to exclude regions inside domain from the pressure 
%           calculation, e.g. solid bodies
%
% blankingfolder: directory containing files with coordinates of the
%                 interface between the blanked region and the flow
%
% blankingroot: root string of filename with coordinates of blanked region
%               boundary (i.e. filename string that precedes file number)
%
% first: first velocity field file number
%
% last: last velocity field file number
%
% increment: integer increment between velocity field numbers to be
%            analyzed
%
% deltaT: time in secconds corresponding to the integer increment set above
%
% numformat: file number format; e.g. inroot_001 is '%03d' format, 
%            blankingroot_00001 is '%05d' format 
%
% fileextension: file extension
%
% separator: character separating columns in data files; use '\t' to 
%            specify tab separator, use '' to specify single space
%            separator
%
% numheaderlines: number of header lines in velocity field data files
%
% nodecrop: number of velocity nodes at domain boundary to blank, e.g. to 
%           remove errant boundary data; value must be >= 1
% 
% lengthcalib_axis: meters per unit of length for velocity field axes 
%                   (e.g. 1/1000 if velocity field axes are in mm)
% 
% lengthcalib_vel: meters per unit of length for velocity field vectors 
%                  (e.g. 1/1000 if velocity field vectors are in mm/sec)
% 
% timecalib_vel: seconds per unit of time for velocity field vectors 
%                (e.g. 60 if velocity field vectors are in mm/min)
% 
% nu: fluid kinematic viscosity in m^2/s
% 
% rho: fluid density in kg/m^3
% 
% gradientplane: for 3D velocity field data, set to 'xy' or 'yz' if 
%                velocity gradients are significantly greater in the 
%                indicated plane; else set to 'iso'. If velocity gradients
%                are significantly greater in the 'xz' plane, relabel the
%                axes so that the plane is 'xy' or 'yz'.
% 
% viscous: set to 1 to include viscous term in pressure gradient 
%          calculation; set to 0 to omit viscous term in pressure gradient 
%          calculation
% 
% smooth_t: set to 1 to implement temporal spline smoothing filter of input 
%           velocity fields
%
% smooth_dp: set to 1 to implement nearest-neighbor smoothing of pressure 
%            gradient fields
% 
% smooth_p: set to 1 to implement nearest-neighbor smoothing of pressure 
%           field
%
% plot_delay: time delay (in seconds) between the display of successive 
%             velocity and pressure fields for time series data sets
%
% plotfigure: set to 1 to plot pressure contours for each velocity field
% 
% export_data: set to 1 to save pressure data to file 
%
% export_format: set to 'ascii' to save an unformatted ASCII data file (one 
%                file per velocity field); set to 'tecplot' to save a 
%                Tecplot-compatible data file (one combined file).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    args = struct(...
      'datafolder'           , './examples/jellyfish', ...
      'inroot'               , 'jellyfish_', ...
      'outroot'              , 'queen2_jellyfish_', ...
      'blanking'             , 0, ...
      'blankingfolder'       , './', ...
      'blankingroot'         , 'interface_', ...
      'first'                , 1, ...
      'last'                 , 100, ...
      'increment'            , 1, ...
      'deltaT'               , 0.005, ...
      'numformat'            , '%05d', ...
      'fileextension'        , '.dat', ...
      'separator'            , ',', ...
      'numheaderlines'       , 0, ...
      'nodecrop'             , 1, ...
      'lengthcalib_axis'     , 1/1000, ...
      'lengthcalib_vel'      , 1, ...
      'timecalib_vel'        , 1, ...
      'nu'                   , 1e-6, ...
      'rho'                  , 1000, ...
      'gradientplane'        , 'yz', ... 
      'viscous'              , 1, ...
      'smooth_t'             , 1, ...
      'smooth_dp'            , 1, ...
      'smooth_p'             , 1, ...
      'plot_delay'           , 1, ...
      'export_data'          , 1, ...
      'export_format'        , 'ascii' ... 
    );