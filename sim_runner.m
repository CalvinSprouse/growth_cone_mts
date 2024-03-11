% use this script to actually run axonal_mts.m
clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHANGE THE PARAMTER IN ALL LOCATIONS BEFORE RUNNING %
% current paramter: lave
%
% checklist:
% (l.37) change_variable   |X|
% (l.45) export_dir        |X|
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define sims per variable
sim_count = 10;

% define variable values
var_val = round(linspace(5, 50, 20));

% define progress bar items
prog_tot = sim_count * length(var_val);
prog_val = 0;
prog_bdr = "____________________________";

% iterate over parameters
for run_i = 1:length(var_val)
	% iterate over sim_count
	for sim_j = 1:sim_count
		% update progress bar
		clc;
		fprintf("%s\nProgress 'Bar' (%.2f%%)\n", prog_bdr, prog_val*100)
		fprintf("i=%d \t j=%d \n%s\n", run_i, sim_j, prog_bdr);
		
		% load defaults
		sim_defaults;

		% change variables
		lave = var_val(run_i);
		ttot = 500;

		% define save locations
		save_str = sprintf("sweep_p%i-s%i.mat", run_i, sim_j);
		export_dir = fullfile("SimExport", "final_lave");
		export_file = fullfile(export_dir, save_str);

		% run simulation
		axonal_mts;

		% update progress tracker		
		prog_val = prog_val + 1/prog_tot;
	end
end

% print final bar (no need, just pretty)
clc;
fprintf("%s\nProgress 'Bar' (%.2f%%)\n", prog_bdr, prog_val*100)
fprintf("i=%d \t j=%d \n%s\n", run_i, sim_j, prog_bdr);
