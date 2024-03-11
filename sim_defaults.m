% sim_defaults

% this file contains the default simulation parameter definitions for use in axonal_mts
% this file will be called automatically by axonal_mts before the simulation starts

% base units:
% time: s,
% length: micron,
% force: pN.

% input parameters constrained by experimental measurements:
% see Craig et. al., MBoC 2017 for more info on these parameters
% protein friction per unit length [pN / um^2]
xi = 0.00144;
% maximum motor density [1/um]
Dmax = 125.0;
% dynein detachment rate per site under zero load [1/s]
doff = 0.37;
% dynein detachment force [pN]
Fdd = 1.74;
% dynein stall force [pN]
Fsd = 2.5;
% dynein forward velocity [um/s]
vfd = 3.5;
% dynein backward velocity [um/s]
vbd = 0.001;

% mt dynamics parameters:
% mt polymerization rate [um/s]
vp = 0.1;
% MT depolymerization rate [um/s]
vm = 0.16;
% catastrophe frequency [1/s]
fpm = 0.0;
% rescue frequency [1/s]
fmp = 0.0;
% rate of switching from stable to dynamic [1/s]
rdyn = 0.0;

% input parameters that are not experimentally constrained:
% axon length [um]
laxon = 100.0;
% mean mt length [um]
lave = 10.0;
% std mt length [um]
lsig = 5;
% dynein attachment rate [1/s]
don = 0.0001;
% xpar attachment rate [1/s]
xonpar = 0.0;
% xanti attachment rate [1/s]
xonantipar = 0.0;
% xpar detachment rate [1/s]
xoffpar = 1.0;
% xanti detachment rate [1/s]
xoffanti = 1.0;
% xpar effective drag coefficient [pN / um]
gammapar = 1000;
% xanti effective drag coefficient [pN / um]
gammaanti = 1000;
% initial mt count [#]
init_nmt = 1;
% initial percent meo [%]
pflipi = 0.25;
% rate of flipping events [1/s]
rflip = 0.0;
% rate of severing events per binding site [1/s]
rsev = 0.0;
% nucleation rate [1/s]
rnuc = 0.0;

% time per step [s]
dt = 0.01;
% total runtime [s]
ttot = 100;

% cell body clearing percentage (%)
% 0 means the mt is cleared on contact
% 1 means the mt is cleared on immersion
% fractional values produce fractional behavior
cb_coeff = 0.5;

% tracks if parameters were already loaded
vars_loaded = 1;

% define the export directory and filename
save_sim = 1;
export_dir = fullfile("SimExport");
export_file = fullfile(export_dir, "axonal_mts.mat");
