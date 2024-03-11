% growth_cone_mts.m

%{
Script description:
This script is a modification to 'Traction_model_052813_3.m' to include
microtubule (mt) action in the growth cone (gc). See 'Craig et al. 2015'
for more information on original model features.

Model description:

%}

% init workspace
clear;

% define model
% 1 means binding occurs at MT tip mostly
% 2 means binding occurs throughout MT
xlink_model = 2;

% configure simulation parameters
% define the system extent or length [micrometer]
L = 12;
% define the number of grid points [ ]
N = 40;
% define the time step [s]
Deltat = 0.0001;
% define length of T domain [micrometer]
LT = 5;

% configure dimensionless free protein parameters
% see Craig et al. 2015 page 8
% measure of actin viscosity
lambda = 3;
% measure of myosin strength
F = 4;
% measure of myoson rate of detachment
gamma = 1000;
% measure of the myoson influence on the adhesion strengthening
b = 15;
% 
Vstar = 0.3;


% calculate grid spacing [micrometer]
Deltax = L / (N-1);
% calculate grid sizes [micrometer]
x = Deltax * (0:N-1);

% calculate yet unkown 'h' terms
% dtdx
dtdx = Deltat / Deltax;
% dtdx^2
dtdx2 = Deltat / (Deltax^2);

% calculate the number of steps
nstep = 100/Deltat;


% calculate protein properties
% veloity [micrometer/s]
v = (N - (1:N)) / (N-1);
v_new = (N - (1:N)) / (N-1);
% traction force [pN]
T = N - (1:N);

% initialize the myosin contractile force array
% F_myo(x) [pN]
Fmyo = ones(size(1:N));

% initialize the myson spacial distribution array
M = ones(1, N);
M_new = ones(1, N);


% loop to reach steady state
for j=1:nstep
    v_new(2:N-1) = v(2:N-1) + dtdx2*(v(3:N) + v(1:N-2) -2*v(2:N-1))+(Deltat/(lambda))*(-T(2:N-1)+Fmyo(2:N-1));

    v_new(N) = v(N-1);
    v = v_new;

    M_new(2:(N)) = M(2:(N))+Deltat*(1-M(2:(N)))-dtdx*gamma*(v(2:(N)).*M(2:(N))-v(1:(N-1)).*M(1:(N-1)));
    M_new(1) = 0;
    M = M_new;

    % update myosin stress
    delm(2:N) = (M(2:(N))-M(1:(N-1)))/Deltax;
    delm(1) = (M(2)-M(1))/Deltax;
    Fmyo = F*M;

    Nadh = (1+b*Fmyo).*exp(-v/Vstar);
    T = Nadh.*v;
end

% mt section
% first we determine the characteristic time scale
% this will be how long it takes actin to cross the network
% characteristic time to reach stability [s]
t_char = Deltax * sum(v.^(-1));

% define some parameters
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
fpm = 0.01;
% rescue frequency [1/s]
fmp = 0.03;
% rate of switching from stable to dynamic [1/s]
rdyn = 0.001;

% input parameters not experimentally constrained:
% dynein attachment rate [1/s]
don = 0.0001;
% xpar attachment rate [1/s]
xonpar = 0.0001;
% xpar detachment rate [1/s]
xoffpar = 0.000;
% xpar effective drag coefficient [pN / um]
gammapar = 1000;
% xanti effective drag coefficient [pN / um]
gammaanti = 1000;
% initial mt count [#]
nmt = 100;

% time per step [s]
dt = 0.01;
% total runtime [s]
ttot = round(t_char);
% calculate steps
nsteps_mts = ttot/dt;

% initialize the mt array (steps, variables, mts)
MTs = zeros(nsteps_mts, 14, nmt);
% 1. initial time (t) [s]
MTs(1, 1, :) = 0.0;
% 2. uniformly randomized initial position for each MT [um]
% MTs(1, 2, :) = L .* rand(1,1,nmt);
MTs(1, 2, :) = 6;
% 3. initial velocity (v) [um / s]
MTs(1, 3, :) = 0.0;
% 4. initial force (Fdp) [pN]
MTs(1, 4, :) = 0.0;
% 5. initial force (Fdm) [pN]
MTs(1, 5, :) = 0.0;
% 6. initial forward-pulling dynein attachment number (Ndp) [#]
MTs(1, 6, :) = 0.0;
% 7. initial backward-pulling dynein attachment number (Ndm) [#]
MTs(1, 7, :) = 0.0;
% 8. initial xpar attachment number (Nxpar) [#]
MTs(1, 8, :) = 0.0;
% 9. initial xanti attachment number (Nxanti) [#]
MTs(1, 9, :) = 0.0;
% 10. initial length = 30 micrometers
MTs(1, 10, :) = 30;
% 11. polarity; 0 is peo, 1 is meo.
MTs(1, 11, :) = 0;
% 12. stability; 0 is stable, 1 is growing, 2 is shrinking
MTs(1, 12, :) = 0;
% 13. force from myosin binding
MTs(1, 13, :) = 0.0;
% 14. The status of MTs as "alive" 1 or "not" 0
MTs(1, 14, :) = 1;

% next step primary loop
for si = 1:(nsteps_mts-1)
	% propogate forward all values from the previous timestep
	MTs(si + 1, :, :) = MTs(si, :, :);

	% update time
	MTs(si + 1, 1, :) = MTs(si, 1, :) + dt;

	% update position
	MTs(si + 1, 2, :) = MTs(si, 2, :) + MTs(si, 3, :)*dt;

	% placeholder arrays for attachment numbers
	ndp = MTs(si, 6, :);
	ndm = MTs(si, 7, :);
	nxpar = MTs(si, 8, :);
	nxanti = MTs(si, 9, :);

	% apply per-mt updates
	for mt_i = 1:nmt
		% update the length of dynamic MTs
		if MTs(si, 12, mt_i) == 1
			% update growing mt length
			MTs(si + 1, 10, mt_i) = MTs(si, 10, mt_i) + vp*dt;
		elseif (MTs(si, 12, mt_i) == 2)
			% update shrinking mt length
			MTs(si + 1, 10, mt_i) = MTs(si, 10, mt_i) - vm*dt;
		end

		% mts may switch between stable and dynamic
		if (rand() < rdyn*dt) && (MTs(si, 12, mt_i) == 0)
			if (rand() < 0.5)
				% stable to growing
				MTs(si + 1, 12, mt_i) = 1;
			else
				% stable to shrinking
				MTs(si + 1, 12, mt_i) = 0;
			end
		elseif (rand() < fpm*dt && MTs(si, 12, mt_i) == 1)
			% catastrophe: growing to shrinking
			MTs(si + 1, 12, mt_i) = 2;
		elseif (rand() < fmp*dt && MTs(si, 12, mt_i) == 2)
			% rescue: shrinking to growing
			MTs(si + 1, 12, mt_i) = 1;
		end

		% define some abbreviated variable names for mt_i at step_i (si)
		Fdp = MTs(si, 4, mt_i);
		Fdm = MTs(si, 5, mt_i);

		% apply detachment events; see Crag et. al., MBoC 2017

		% detachment rate for forward-pulling dynein motors
		if ndp(mt_i) == 0
			% no motors
			roffdp = 0;
		elseif (Fdp < ndp(mt_i)*Fsd)
			% load < stall force
			roffdp = doff*exp(Fdp./(ndp(mt_i)*Fdd))*dt;
		else
			% load >= stall force
			roffdp = doff*(4/(1-exp(-Fdp./(ndp(mt_i)*1.97))))*dt;
		end

		% detachment rate for backward-pulling dynein motors
		if ndm(mt_i) == 0
			roffdm = 0;
		elseif (Fdm < ndm(mt_i)*Fsd)
			% load < stall force
			roffdm = doff*exp(Fdm./(ndm(mt_i)*Fdd))*dt;
		else
			% load > stall force
			roffdm = doff*(4/(1-exp(-Fdm./(ndm(mt_i)*1.97))))*dt;
		end

		% detachment rate for xpar
		roffxpar = xoffpar*dt;

		% allow forward-pulling dynein detachment
		ndp(mt_i) = ndp(mt_i) - nnz(rand(1, ndp(mt_i)) < roffdp);

		% ensure forward-pulling dynein number does not go negative
		if ndp(mt_i) < 0
			ndp(mt_i) = 0;
		end

		% allow backward-pulling dynein detachment
		ndm(mt_i) = ndm(mt_i) - nnz(rand(1, ndm(mt_i)) < roffdm);

		% ensure backward-pulling dynein number does not go negative
		if ndm(mt_i) < 0
			ndm(mt_i) = 0;
		end

		% allow xpar detachment
		nxpar(mt_i) = nxpar(mt_i) - nnz(rand(1, nxpar(mt_i)) < roffxpar);

		% ensure xpar number does not go negative
		if nxpar(mt_i) < 0
			nxpar(mt_i) = 0;
		end

		% ensure no xpar remain if outside L
		if MTs(si+1, 2, mt_i) > L
			nxpar(mt_i) = 0;
		end

		% calculate attachment rates
		% begin with the dynein zone
		% calculate overlap with C domain
		overlap = ( MTs(si,2,mt_i) + MTs(si,10,mt_i) ) - ( L+LT );

		% ensure overlap is non-negative
		if overlap < 0
			overlap = 0;
		end

		% calculate the total number of binding sites
		Ntoti = round(Dmax*MTs(si, 10, mt_i));

		% sum the motors on the mt
		mt_i_motors = MTs(si, 6, mt_i) + MTs(si, 7, mt_i) + MTs(si, 8, mt_i) + MTs(si, 9, mt_i);

		% estimate the number of overlapping binding sites
		Noverlap = (overlap/MTs(si, 10, mt_i))*(Ntoti - mt_i_motors);
		if Noverlap < 0
			Noverlap = 0;
		end

		% on-rate for dynein
		rond = Noverlap * don * dt;

		% now do the same for cross linkers to the actomyosin
		% calculate overlap with the P domain
		overlap = L - MTs(si,2,mt_i);

		% ensure overlap is non-negative
		if overlap < 0
			overlap = 0;
		end

		% calculate the total number of binding sites
		Ntoti = round(Dmax*MTs(si, 10, mt_i));

		% sum the motors on the mt
		mt_i_motors = MTs(si, 6, mt_i) + MTs(si, 7, mt_i) + MTs(si, 8, mt_i) + MTs(si, 9, mt_i);

		% estimate the number of overlapping binding sites
		Noverlap = (overlap/MTs(si, 10, mt_i))*(Ntoti - mt_i_motors);
		if Noverlap < 0
			Noverlap = 0;
		end

		% on rate for xlinker
		ronxpar = Noverlap * xonpar * dt;

		% permit attachments for dynein
		if rand() < rond
			if rand() < 0.5
				% forward dynein
				ndp(mt_i) = ndp(mt_i) + 1;
			else
				%backward dynein
				ndm(mt_i) = ndm(mt_i) + 1;
			end
		end

		% permit xpar attach
		if rand() < ronxpar
			nxpar(mt_i) = nxpar(mt_i) + 1;
		end

		% assign updated attachment numbers for mt i
		MTs(si + 1, 6, mt_i) = ndp(mt_i);
		MTs(si + 1, 7, mt_i) = ndm(mt_i);
		MTs(si + 1, 8, mt_i) = nxpar(mt_i);
		MTs(si + 1, 9, mt_i) = nxanti(mt_i);

		% sum the new counts of attached proteins
		attachment_number = ndp(mt_i) + ndm(mt_i) + nxpar(mt_i) + nxanti(mt_i);

		% find the forces from the new xlinkers
		if xlink_model == 1
			% xlinkers bounded at the tip
			approx_n = round((abs(MTs(si, 2, mt_i))*N/L) + 1);
			if approx_n > 40
				approx_n = 40;
			end
			MTs(si + 1, 13, mt_i) = nxpar(mt_i) * Fmyo(approx_n);
			cross_v = v(approx_n);
		elseif xlink_model == 2
			% xlinkers bounded evenly
			approx_n_min = round((abs(MTs(si, 2, mt_i))*N/L) + 1);
			n_range = round( approx_n_min:nxpar(mt_i):N );
			MTs(si+1, 13, mt_i) = sum(Fmyo(n_range));
			cross_v = sum(v(n_range))/nxpar(mt_i);
			if nxpar(mt_i) == 0
				cross_v = 0;
			end
		end

		% calculate the effective number of minus motors
		np = ndp(mt_i) + nxpar(mt_i);

		% define forward and backward velocity
		% use weighted averages
		vb = vbd;
		vf = (ndp(mt_i)*vfd + nxpar(mt_i)*cross_v)/np;
		if np == 0
			vf = vfd;
		end

		% update v, Fdp, for time i+1 based on new attachment numbers
		% plus pulling dominant elseif minus pulling dominant
		if np > ndm(mt_i)
			% calculate v
			% v_unsigned = vfd*(ndp(mt_i) - ndm(mt_i))./(ndp(mt_i) + (vfd/vbd)*ndm(mt_i) + (vfd/Fsd)*(xi + nxpar(mt_i)*gammapar + nxanti(mt_i)*gammaanti));
			v_unsigned = vf*(np-ndm(mt_i))./(np+(vf/vb)*ndm(mt_i)+(vf/Fsd)*(xi+nxpar(mt_i)*gammapar+nxanti(mt_i)*gammaanti));

			% if plus-oriented
			if (MTs(si + 1, 11, mt_i) == 0)
				% calculate new velocity
				MTs(si + 1, 3, mt_i) = v_unsigned;
			else
				MTs(si + 1, 3, mt_i) = -v_unsigned;
			end

			% update Fdp
			% first calculate numerator
			Fnum = (np.*(ndm(mt_i)*Fsd*(1 + (vf/vb)) + vf*np.*(xi + nxpar(mt_i)*gammapar + nxanti(mt_i)*gammaanti)));

			% then calculate denominator
			Fdenom = ndp(mt_i) + (vf/vb)*ndm(mt_i) + (vf/Fsd)*(xi + nxpar(mt_i)*gammapar + nxanti(mt_i)*gammaanti);

			% full expression for updated Fdp
			MTs(si + 1, 4, mt_i) = Fnum/Fdenom;

			% update Fdm
			MTs(si + 1, 5, mt_i) = MTs(si + 1, 4, mt_i) - ((xi + nxpar(mt_i)*gammapar + nxanti(mt_i)*gammaanti)*MTs(si + 1, 3, mt_i));
		elseif ndm(mt_i) > np
			% calculate v
			% v_unsigned = -vfd*(ndm(mt_i) - ndp(mt_i))./(ndm(mt_i) + (vfd/vbd)*ndp(mt_i) + (vfd/Fsd)*(xi + nxpar(mt_i)*gammapar + nxanti(mt_i)*gammaanti));
			v_unsigned = -vf*(ndm(mt_i)-np)./(ndm(mt_i)+(vf/vb)*np+(vf/Fsd)*(xi+nxpar(mt_i)*gammapar+nxanti(mt_i)*gammaanti));

			% if MT i has a plus polarity
			if (MTs(si + 1, 11, mt_i) == 0)
				MTs(si + 1, 3, mt_i) = v_unsigned;
			else
				MTs(si + 1, 3, mt_i) = -v_unsigned;
			end

			% update Fdp
			% first calculate numerator
			Fnum = (1 + (vf/vb))*np*ndm(mt_i)*Fsd + vf*ndm(mt_i)*(xi + nxpar(mt_i)*gammapar + nxanti(mt_i)*gammaanti);

			% then calculate denominator
			Fdenom = ndm(mt_i) + (vf/vb)*np + (vf/Fsd)*(xi + nxpar(mt_i)*gammapar + nxanti(mt_i)*gammaanti);

			% full expression for updated Fdm because backwards
			MTs(si + 1, 5, mt_i) = Fnum/Fdenom;

			% update Fdp
			MTs(si + 1, 4, mt_i) = MTs(si + 1, 5, mt_i) + ((xi + nxpar(mt_i)*gammapar + nxanti(mt_i)*gammaanti)*MTs(si + 1, 3, mt_i));

			% no need to check for MT pressing into the distill end since
			% its moving in the wrong direction to do so
		else
			% no force so no velocity
			MTs(si + 1, 3, mt_i) = 0;
			MTs(si + 1, 4, mt_i) = 0;
			MTs(si + 1, 5, mt_i) = 0;
		end
	end
end


% find the index corresponding to maximum traction
[maxT, maxi] = max(T);

% find the corresponding position and velocity to max traction
xmax = x(maxi);
vmax = v(maxi);

% calculate output quantities
% x = Deltax*(0:(N-1));
Tv = v .* exp(-v / Vstar);
Tf = 1 + b*Fmyo;

% velocity, myosin force, and traction spacial distribution figure
% figure
% plot(x, v, 'm', x, Fmyo, 'w--', x, T, 'r', LineWidth=2);
% legend('Actin flow speed','Myosin density','Traction');

% save simulation output to .mat file
A = horzcat(x', v', T', M', Tv', Tf', Nadh');
Parms = vertcat(lambda, F, gamma, b, Vstar, xmax, vmax);
% Parms = vertcat(v0,mu0,F0,T0,Vstar,moff,beta,xmax,vmax);
save("growth_cone_mts_m2", "-v7", "-nocompression");
