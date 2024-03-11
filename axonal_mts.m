% axonal_mts

% brief model description:
% Axonal microtubules slide in response to motor-bases forces from dynein.
% When Dynein cross-links to two MTs, one MT is attached to cargo domain
% and other is attached to motor domains, and the two MTs are slid in opposite
% directions. Non-motile cross-linkers can also crosslink two MTs and induce
% viscous-like resistance relative to sliding. Short MTs that are not attached
% to any other MTs have some chance to flip. All MTs start as stable, not growing
% or shrinking, but can switch into states of growing and shrinking known as
% dynamic instability. See Craig et. al., MBoC 2017 for more information.
%
% this version of the code is modified from "axonal_mts_080822.m" to include:
% a lot.
%
% accroyms/aabbreviations:
% mt = microtubule
% peo = plus-end-out
% meo = minus-end-out
% xpar = parallel cross-linker
% xanti = anti-parallel cross-linker
% lhs = left hand side
% rhs = right hand side

% load the parameters from a parameter definition script unless vars were already loaded
if ~exist("vars_loaded", "var")
	clear;
	sim_params;
end

% calculate number of steps [#]
nstep = ttot/dt;

% initialize a polarity array for nucleation events
polarity = zeros(laxon + 1, 3);

% track the total number of MTs
nmt = init_nmt;

% generate initial mt polaritys using a randomly sorted logical array
init_pol = zeros([1, 1, nmt]);
init_pol(1 : round(nmt * pflipi)) = 1;
init_pol = init_pol(randperm(length(init_pol)));

% generate initial mt positions using a normal distribution about the axon
init_pos = laxon * rand(1,1,nmt);

% generate initial mt lengths
% TODO: handle negative lengths
init_len = normrnd(lave, lsig, 1, 1, nmt);

% initialize the mt array (steps, variables, mts)
MTs = zeros(nstep, 14, nmt);
% 1. initial time (t) [s]
MTs(1, 1, :) = 0.0;
% 2. uniformly randomized initial position for each MT [um]
MTs(1, 2, :) = init_pos;
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
% 10. randomized initial length
MTs(1, 10, :) = 2*lave*rand([1, 1, nmt]);
% 11. polarity; 0 is peo, 1 is meo.
MTs(1, 11, :) = init_pol;
% 12. stability; 0 is stable, 1 is growing, 2 is shrinking
MTs(1, 12, :) = 0;
% 13. initial terminal force (Fte) [pN]
MTs(1, 13, :) = 0.0;
% 14. The status of MTs as "alive" 1 or "not" 0
MTs(1, 14, :) = 1;


% program primary loop
for si = 1:(nstep-1)
	% propogate forward all values from the previous timestep
	MTs(si + 1, :, :) = MTs(si, :, :);

	% update time
	MTs(si + 1, 1, :) = MTs(si, 1, :) + dt;

	% update position
	MTs(si + 1, 2, :) = MTs(si, 2, :) + MTs(si, 3, :)*dt;

	% calculate polarity pattern for nucleation events
	polarity = calc_polarity_pattern(MTs, si, laxon, nmt);

	% placeholder arrays for attachment numbers
	ndp = MTs(si, 6, :);
	ndm = MTs(si, 7, :);
	nxpar = MTs(si, 8, :);
	nxanti = MTs(si, 9, :);

	% apply per-mt updates and calculations (slow part)
	for mt_i = 1:nmt
		% check if MT is alive, if not then skip
		if ~MTs(si, 14, mt_i)
			continue;
		end

		% check if MT is some percentage within the cell body
		if (MTs(si, 2, mt_i) + cb_coeff*MTs(si + 1, 10, mt_i)) < 0
			% "kill" the mt
			MTs(si + 1, 14, mt_i) = 0;
		end

		% update the length of dynamic MTs
		if MTs(si, 12, mt_i) == 1
			% update growing mt length
			MTs(si + 1, 10, mt_i) = MTs(si, 10, mt_i) + vp*dt;
		elseif (MTs(si, 12, mt_i) == 2)
			% update shrinking mt length
			MTs(si + 1, 10, mt_i) = MTs(si, 10, mt_i) - vm*dt;
		end

		% if mt length shrinks to 0
		if MTs(si + 1, 10, mt_i) <= 0
			% "kill" current MT
			MTs(si + 1, 14, mt_i) = 0;
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
		t = MTs(si, 1, mt_i);
		x = MTs(si, 2, mt_i);
		v = MTs(si, 3, mt_i);
		Fdp = MTs(si, 4, mt_i);
		Fdm = MTs(si, 5, mt_i);
		L = MTs(si, 10, mt_i);
		P = MTs(si, 11, mt_i);

		% apply detachment events; see Crag et. al., MBoC 2017

		% detachment rate for forward-pulling dynein motors
		if ndp(mt_i) == 0
			% no motors
			roddfp = 0;
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

		% detachment rate for xanti
		roffxanti = xoffanti*dt;

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

		% allow xanti detachment
		nxanti(mt_i) = nxanti(mt_i) - nnz(rand(1, nxanti(mt_i)) < roffxanti);

		% ensure xanti number does not go negative
		if nxanti(mt_i) < 0
			nxanti(mt_i) = 0;
		end

		% calculate attachment rates
		for mt_k = (mt_i + 1):nmt
			% check that comparison mt is alive
			if ~MTs(si + 1, 14, mt_k)
				continue
			end

			% find the overlap between mt_i and mt_k
			if MTs(si, 2, mt_k) > MTs(si, 2, mt_i)
				% lhs of mt_k > lhs of mt_i
				overlap = MTs(si, 2, mt_i) + MTs(si, 10, mt_i) - MTs(si, 2, mt_k);
			else
				% rhs of mt k < lhs of mt i
				overlap = MTs(si, 2, mt_k) + MTs(si, 10, mt_k) - MTs(si, 2, mt_i);
			end

			% ensure overlap is non-negative
			if overlap < 0
				overlap = 0;
			end

			% calculate the total number of binding sites
			Ntoti = round(Dmax*MTs(si, 10, mt_i));
			Ntotk = round(Dmax*MTs(si, 10, mt_k));

			% sum the motors on each MT
			mt_i_motors = MTs(si, 6, mt_i) + MTs(si, 7, mt_i) + MTs(si, 8, mt_i) + MTs(si, 9, mt_i);
			mt_k_motors = MTs(si, 6, mt_k) + MTs(si, 7, mt_k) + MTs(si, 8, mt_k) + MTs(si, 9, mt_k);

			% estimate the number of overlapping binding sites for i and k
			Noverlapi = (overlap/MTs(si, 10, mt_i))*(Ntoti - mt_i_motors);
			Noverlapk = (overlap/MTs(si, 10, mt_k))*(Ntotk - mt_k_motors);

			% take the smaller of the two to calculate attachment rates
			Noverlap = min([Noverlapk, Noverlapi]);
			if (Noverlap < 0)
				Noverlap = 0;
			end

			% calculate on-rate for each population for i and k
			% xpar and xanti are exclusive by relative polarity
			if (MTs(si, 11, mt_i) == MTs(si, 11, mt_k))
				ronxpar = Noverlap * xonpar * dt;
				ronxanti = 0;
			else
				ronxpar = 0;
				ronxanti = Noverlap * xonantipar * dt;
			end

			% on-rate for dynein
			rond = Noverlap * don * dt;

			% allow motors and cross-linkers to attach to mt j and k
			if rand() < rond
				% if parallel then one gets forward and one gets backward
				if MTs(si, 11, mt_i) == MTs(si, 11, mt_k)
					if (rand() < 0.5)
						% mt_i gains forward mt_k gains backwards
						ndp(mt_i) = ndp(mt_i) + 1;
						ndm(mt_k) = ndm(mt_k) + 1;
					else
						% mt_k gains forward mt_i gains backwards
						ndp(mt_k) = ndp(mt_k) + 1;
						ndm(mt_i) = ndm(mt_i) + 1;
					end
				else
					% if anti-parallel then both gain a forward motor
					ndp(mt_i) = ndp(mt_i) + 1;
					ndp(mt_k) = ndp(mt_k) + 1;
				end
			end

			% both mts gain xpar
			if (rand < ronxpar)
				nxpar(mt_i) = nxpar(mt_i) + 1;
				nxanti(mt_k) = nxanti(mt_k) + 1;
			end

			% both mts gain xanti
			if (rand < ronxanti)
				nxanti(mt_i) = nxanti(mt_i) + 1;
				nxanti(mt_k) = nxanti(mt_k) + 1;
			end
		end

		% assign updated attachment numbers for mt i
		% mt k will have updated attachment numbers on its own for loop
		MTs(si + 1, 6, mt_i) = ndp(mt_i);
		MTs(si + 1, 7, mt_i) = ndm(mt_i);
		MTs(si + 1, 8, mt_i) = nxpar(mt_i);
		MTs(si + 1, 9, mt_i) = nxanti(mt_i);

		% sum the new counts of attached proteins
		attachment_number = ndp(mt_i) + ndm(mt_i) + nxpar(mt_i) + nxanti(mt_i);

		% allow short mts to flip orientation
		% short means < 2 microns length and no motor attachments
		% TODO: Replace with cylindrical diffusion rate
		if (MTs(si + 1, 10, mt_i) < 2.0 && attachment_number == 0)
			if (rand < rflip*dt)
				% plus to minus else minus to plus
				if (MTs(si + 1, 11, mt_i) == 0)
					MTs(si + 1, 11, mt_i) = 1;
				else
					MTs(si + 1, 11, mt_i) = 0;
				end
			end
		end

		% update v, Fdp, for time i+1 based on new attachment numbers
		% plus pulling dominant elseif minus pulling dominant
		if ndp(mt_i) > ndm(mt_i)
			% calculate v
			v_unsigned = vfd*(ndp(mt_i) - ndm(mt_i))./(ndp(mt_i) + (vfd/vbd)*ndm(mt_i) + (vfd/Fsd)*(xi + nxpar(mt_i)*gammapar + nxanti(mt_i)*gammaanti));

			% if plus-oriented
			if (MTs(si + 1, 11, mt_i) == 0)
				% calculate new velocity
				MTs(si + 1, 3, mt_i) = v_unsigned;
			else
				MTs(si + 1, 3, mt_i) = -v_unsigned;
			end

			% update Fdp
			% first calculate numerator
			Fnum = (ndp(mt_i).*(ndm(mt_i)*Fsd*(1 + (vfd/vbd)) + vfd*ndp(mt_i).*(xi + nxpar(mt_i)*gammapar + nxanti(mt_i)*gammaanti)));

			% then calculate denominator
			Fdenom = ndp(mt_i) + (vfd/vbd)*ndm(mt_i) + (vfd/Fsd)*(xi + nxpar(mt_i)*gammapar + nxanti(mt_i)*gammaanti);

			% full expression for updated Fdp
			MTs(si + 1, 4, mt_i) = Fnum/Fdenom;

			% update Fdm
			MTs(si + 1, 5, mt_i) = MTs(si + 1, 4, mt_i) - ((xi + nxpar(mt_i)*gammapar + nxanti(mt_i)*gammaanti)*MTs(si + 1, 3, mt_i));
		elseif (ndm(mt_i) < ndp(mt_i))
			% calculate v
			v_unsigned = -vfd*(ndm(mt_i) - ndp(mt_i))./(ndm(mt_i) + (vfd/vbd)*ndp(mt_i) + (vfd/Fsd)*(xi + nxpar(mt_i)*gammapar + nxanti(mt_i)*gammaanti));

			% if MT i has a plus polarity
			if (MTs(si + 1, 11, mt_i) == 0)
				MTs(si + 1, 3, mt_i) = v_unsigned;
			else
				MTs(si + 1, 3, mt_i) = -v_unsigned;
			end

			% update Fdp
			% first calculate numerator
			Fnum = (1 + (vfd/vbd))*ndp(mt_i)*ndm(mt_i)*Fsd + vfd*ndm(mt_i)*(xi + nxpar(mt_i)*gammapar + nxanti(mt_i)*gammaanti);

			% then calculate denominator
			Fdenom = ndm(mt_i) + (vfd/vbd)*ndp(mt_i) + (vfd/Fsd)*(xi + nxpar(mt_i)*gammapar + nxanti(mt_i)*gammaanti);

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

		% check if MT is pressing into the distill end / growth cone
		% that is if the MT is touching the distill end and has positive velocity
		if ((MTs(si + 1, 2, mt_i) + MTs(si + 1, 10, mt_i) > laxon) && (MTs(si + 1, 3, mt_i) > 0))
			% then the MT will have no velocity
			MTs(si + 1, 3, mt_i) = 0;

			% regardless of polarity it seems that
			% Fdp - Fdm = friction
			% the frictional force is then the force that
			% prevents the MT from moving
			% but since the MT has no velocity that
			% force must come from the hardwall not drag
			% since Fdp is either pos 4 or 5 depending on
			% polarity we just take the absolute value
			% it makes no sense for the wall to exert a negative force
			MTs(si + 1, 13, mt_i) = abs(MTs(si + 1, 4, mt_i) - MTs(si + 1, 5, mt_i));
		end

		% allow severing events to occur
		% first we ask "will a severing event occur on this MT"
		% then we ask "where along the MT will it occur"
		% then the MT is seperated into a "left" and "right" piece
		% the left piece remains as the original MT
		% the right piece becomes a new MT
		% the motor count on each will be the percentage of
		% original MT length that remains
		% if the mt is seperated with 25% on left and 75% right
		% then 25% of the original motors stay on the left
        sever_chance = MTs(si+1, 10, mt_i) * rsever * dt;
		if (rand() < sever_chance)
			% then we ask where along the MT will this occur
			sever_pos = rand()*MTs(si + 1, 10, mt_i);

			% calculate the left and right percentage of MT
			left_percent = sever_pos / MTs(si + 1, 10, mt_i);
			right_percent = 1 - left_percent;

			% reduce the current MT to the left percent
			MTs(si + 1, 10, mt_i) = left_percent*MTs(si + 1, 10, mt_i);

			% create a new MT that duplicates the old mt and increase NMT
			MTs(:, :, nmt + 1) = MTs(:, :, mt_i);

			% ensure the new MT is marked as "dead" at times before current
			MTs(:, 14, nmt + 1) = 0;

			% since values are passed this should be fine, we just make
			% the MT all "dead" then start as alive, future time steps
			% should "pass" this value
			MTs(si + 1, 14, nmt + 1) = 1;

			% update the motor count of the left MT
			% columns 6, 7, 8, 9
			MTs(si + 1, 6, mt_i) = round(MTs(si + 1, 6, mt_i)*left_percent);
			MTs(si + 1, 7, mt_i) = round(MTs(si + 1, 7, mt_i)*left_percent);
			MTs(si + 1, 8, mt_i) = round(MTs(si + 1, 8, mt_i)*left_percent);
			MTs(si + 1, 9, mt_i) = round(MTs(si + 1, 9, mt_i)*left_percent);

			% update the motor count of the right MT
			MTs(si + 1, 6, nmt + 1) = round(MTs(si + 1, 6, nmt + 1)*right_percent);
			MTs(si + 1, 7, nmt + 1) = round(MTs(si + 1, 7, nmt + 1)*right_percent);
			MTs(si + 1, 8, nmt + 1) = round(MTs(si + 1, 8, nmt + 1)*right_percent);
			MTs(si + 1, 9, nmt + 1) = round(MTs(si + 1, 9, nmt + 1)*right_percent);

			% ensure all motor values are > 0 (this is long)
			if (MTs(si + 1, 6, mt_i) < 0)
				MTs(si + 1, 6, mt_i) = 0;
			end

			if (MTs(si + 1, 7, mt_i) < 0)
				MTs(si + 1, 7, mt_i) = 0;
			end

			if (MTs(si + 1, 8, mt_i) < 0)
				MTs(si + 1, 8, mt_i) = 0;
			end

			if (MTs(si + 1, 9, mt_i) < 0)
				MTs(si + 1, 9, mt_i) = 0;
			end

			if (MTs(si + 1, 6, nmt + 1) < 0)
				MTs(si + 1, 6, nmt + 1) = 0;
			end

			if (MTs(si + 1, 7, nmt + 1) < 0)
				MTs(si + 1, 7, nmt + 1) = 0;
			end

			if (MTs(si + 1, 8, nmt + 1) < 0)
				MTs(si + 1, 8, nmt + 1) = 0;
			end

			if (MTs(si + 1, 9, nmt + 1) < 0)
				MTs(si + 1, 9, nmt + 1) = 0;
			end

			% iterate nmt
			nmt = nmt + 1;
		end

		% allow nucleation events to occur
		% nucleation events occur only at other MTs
		% somewhere along their length
		% new MT polarity is weighted random based on
		% nearby MTs and starts small
		if (rand() < (rmtnucleate/nmt)*dt)
			% create a new MT with the parent MT history
			MTs(:, :, nmt + 1) = MTs(:, :, mt_i);

			% mark it as dead at all steps
			MTs(:, 14, nmt + 1) = 0;

			% mark it as alive at current step
			% this should carry over
			MTs(si + 1, 14, nmt + 1) = 1;

			% give it a short length, 0.1 microns
			MTs(si + 1, 10, nmt + 1) = 0.1;

			% place it randomly along the parent tube
			% roll for "where" along the parent tube the mt started nucleating
			MTs(si + 1, 2, nmt + 1) = MTs(si + 1, 2, nmt + 1) + rand()*MTs(si + 1, 10, mt_i);

			% make it grow
			MTs(1, 12, :) = 1;

			% reset velocity
			MTs(si + 1, 3, nmt + 1) = 0.0;

			% reset the protein counts
			MTs(si + 1, 4, nmt + 1) = 0.0;
			MTs(si + 1, 5, nmt + 1) = 0.0;
			MTs(si + 1, 6, nmt + 1) = 0.0;
			MTs(si + 1, 7, nmt + 1) = 0.0;
			MTs(si + 1, 8, nmt + 1) = 0.0;
			MTs(si + 1, 9, nmt + 1) = 0.0;

			% randomize the polarity based on regional
			% polarity pattern
			loc = round(MTs(si + 1, 2, nmt + 1));

			% check if loc outside the axon
			% this is when a distill touching MT nucleates
			% AND rolls nearly the maximum length adjustment
			if (loc > laxon)
				loc = round(laxon);
			end

			% check for negative loc???
			if (loc < 0)
				loc = 1;
			end

			% get, from polarity, the odds to roll against
			fmin = polarity(loc + 1, 3);

			% roll against fmin
			if (rand() < fmin)
				MTs(si + 1, 11, nmt + 1) = 1;
			else
				MTs(si + 1, 11, nmt + 1) = 0;
			end

			% update the nmt counter
			nmt = nmt + 1;
		end

		% allow random nucleation of an mt
		if (rand() < rnuc*dt)
			% create a new mt at a random location along the axon
			MTs(:, :, nmt + 1) = MTs(:, :, nmt);

			% mark it as dead at all steps
			MTs(:, 14, nmt + 1) = 0;

			% mark it as alive at the current step
			MTs(si + 1, 14, nmt + 1) = 1;

			% give it a short length
			MTs(si + 1, 10, nmt + 1) = 0.1;

			% place it randomly along the axon
			MTs(si + 1, 2, nmt + 1) = rand()*laxon;

			% reset velocity
			MTs(si + 1, 3, nmt + 1) = 0.0;

			% reset the protein counts
			MTs(si + 1, 4, nmt + 1) = 0.0;
			MTs(si + 1, 5, nmt + 1) = 0.0;
			MTs(si + 1, 6, nmt + 1) = 0.0;
			MTs(si + 1, 7, nmt + 1) = 0.0;
			MTs(si + 1, 8, nmt + 1) = 0.0;
			MTs(si + 1, 9, nmt + 1) = 0.0;

			% make growing
			MTs(1, 12, :) = 1;

			% randomiz the polarity based on regional pattern
			loc = round(MTs(si + 1, 2, nmt + 1));

			% check if loc outside the axon
			% this is when a distill touching MT nucleates
			% AND rolls nearly the maximum length adjustment
			if (loc > laxon)
				loc = round(laxon);
			end

			% check for negative loc???
			if (loc < 0)
				loc = 1;
			end

			% get, from polarity, the odds to roll against
			fmin = polarity(loc + 1, 3);

			% roll against fmin
			if (rand() < fmin)
				MTs(si + 1, 11, nmt + 1) = 1;
			else
				MTs(si + 1, 11, nmt + 1) = 0;
			end

			% update the nmt counter
			nmt = nmt + 1;
		end
	end
end

% save the simulation
if save_sim
	[~, ~, ~] = mkdir(export_dir);
	save(export_file, "-nocompression", "-v7");
end