% Brief description of model in this program:

% Axonal microtubules sliding in response to motor-bases forces from dynein
% When dynein cross-links two MTs, one MT is attached to cargo domain and
% other is attached to motor domains, and the two MTs are slid in opposite directions. 
%Non-motile cross-linkers can also crosslink two MTs and induce viscous-like resistance to relative slideing. 
%Short MTs that are not attached to any other MT have some chance to flip.
% All MTs start as 'stable': not growing or shrinking, but can switch into
% states of growing and shrinking known as dynamic instability. 

% Input parameters:
% Base unit of time: s, Base unit of length: micron, Base unit of force: pN

% Input parameters constrained by experimental measurements:
% See Craig et al, MBoC 2017 for more info on these parameters
xi = 0.00144; % Protein friction per unit length (pNs/micron^2)
Dmax = 125.0; % Maximum motor density (#/micron)
doff = 0.37; % Dynein detachment rate under zero load (1/s) 
Fdd = 1.74; % Dynein detachment force (pN)
koff= 0.1; % Kinesin detachment rate under zero load (1/s)
Fdk=1.5; % Kinesin detachment force (pN)
Fs = 2.5; % Stall force (pN) (Assumption: Fs=Fsd=Fsk)
vf = 3.5; % Motor forward velocity (microns/s) (Assumption: v=vfd=vfk)
vb = 0.001; % Motor backward velocity (microns/s) (Assumption: vbd=vbk)

% MT dynamics parameters, can set based on experimental measurements
vp=0.1; % MT polymerization speed (micron/s)
vm=0.16; % MT depolymerization speed (micron/s)
fpm=0.01; % catastrophe frequency (1/s)
fmp=0.03; % rescue frequency (1/s)
rdyn=0.0; % rate of switching from stable to dynamic. If set to zero, MT lengths do not change. If set to 1/dt all MTs become dynamic

% Input parameters that are not experimentally constrained
Laxon= 60.0; % Length of axon
Lave = 10.0; % Average microtubule length (microns)
don = 0.001; % Dynein attachment rate (1/s)
kon = 0.000; % Kinesin attachment rate (1/s)
xonpar = 0.0000; % Cross-linker attachment rate (1/s) for parallel MTs
xonantipar = 0.0000; % Cross-linker attachment rate (1/s) for anti-parallel MTs
xoffpar = 0.01; % Cross-linker detachment rate (1/s)
xoffanti = 0.01; % Cross-linker detachment rate (1/s)
gammapar = 1000; % Cross-linker effective drag coefficient (parallel MTs) (pNs/micron)
gammaanti = 1000; % Cross-linker effective drag coefficient (anti-parallel MTs) (pNs/micron)
NMT = 1000; % Number of MTs
pflipi = 0.5; % probability to flip during initialization
rflip = 0.0; %flipping rate per second during main program

% Program time step
dt=0.01;
ttot=10;
ntime=ttot/dt;

% Initialize variables
Numcleared=0;
Numflip=pflipi*NMT;
A=zeros(ntime+1,14,NMT); % Rows are different times; Columns are: t,x,v,Fdp,Fdm,Nd1,Nd2,Nd3,Nk1,Nk2,Nk3,Nxpar,Nxanti,L,P
A(1,1,:) = 0.0; %Initial time (t)
A(1,2,:)=Laxon*rand([1,1,NMT]); % Initial position (x)
%A(1,2,1:Numflip)=Laxon*rand([1,1,Numflip]); % Randomized initial position for each flipped MT
%A(1,2,(Numflip+1):NMT)=Laxon*rand([1,1,(NMT-Numflip)]);
A(1,3,:)=0.0; % Initial velocity (v)
A(1,4,:)=0.0; % Initial force (Fdp)
A(1,5,:)=0.0; % Initial force (Fdm)
A(1,6,:)=0.0; % Initial Nd1 (dynein attached in MT sliding config.)
A(1,7,:)=0.0; % Initial Nd2 (dynein attached in cargo config., walking on parallel MT)
A(1,8,:)=0.0; % Initial Nd3 (dynein attached in cargo config., walking on anti-parallel MT)
A(1,9,:)=0.0; % Initial Nk1 (kinesin attached in MT sliding config.)
A(1,10,:)=0.0; % Initial Nk2 (kinesin attached in cargo config., walking on parallel MT)
A(1,11,:)=0.0; % Initial Nk3 (kinesin attached in cargo config., walking on anti-parallel MT)
A(1,12,:)=0.0; % Initial cross-linker attachment number, parallel MTs (Nxpar)
A(1,13,:)=0.0; % Initial cross-linker attachment number, anti-parallel MTs (Nxanti)
A(1,14,:)=2*Lave*rand([1,1,NMT]); % length
% A(1,10,:) designates polarity. P = 0 is 'plus-end-out', P=1 is 'minus-end-out'
%A(1,15,1:Numflip)=1;
%A(1,15,(Numflip+1):NMT)=0;
for j=1:NMT
   if(rand<pflipi & A(1,14,j)<2) % if random number < prob of initial flipping and MT length < 2micron, then make this MT minus-end-out
        A(1,15,j)=1; % minus end out
    else
        A(1,15,j)=0; % plus end out
    end
end

A(1,16,:)=0; % 0 = stable, 1 = dynamic & growing, 2= dynamic and shrinking


    
    forwardi=A(1,:,A(1,15,:)==0);
    backwardi=A(1,:,A(1,15,:)==1);
    forwardshorti=forwardi(1,:,forwardi(1,14,:)<2);
    backwardshorti=backwardi(1,:,backwardi(1,14,:)<2);
    
    
    figure
    histogram(forwardi(1,2,:))
    hold
    histogram(backwardi(1,2,:))
    xlim([0 Laxon])
    ylim([0 100])
    title('Initial distribution, all MTs')
    ylabel('Total number of MTs')
    xlabel('Distance from cell body, \mum')
    legend('plus-end-out','minus-end-out')
 
  % Calculate percent of + out and - out MTs in each region for bar plot

%Proximal region, defined as 0<x<Laxon/3
nforwardproxi=numel(forwardi(1,2,forwardi(1,2,:)>0 & forwardi(1,2,:)<Laxon/3)); % # plus-out MTs in this region
nbackwardproxi=numel(backwardi(1,2,backwardi(1,2,:)>0 & backwardi(1,2,:)<Laxon/3)); % # of minus-out MTs in this region
percentbackproxi=(nbackwardproxi/(nforwardproxi+nbackwardproxi))*100; % percent with minus-out

% Distal region, defined as x>(2/3)*Laxon
nforwarddisi=numel(forwardi(1,2,forwardi(1,2,:)>(2/3)*Laxon & forwardi(1,2,:)<Laxon));
nbackwarddisi=numel(backwardi(1,2,backwardi(1,2,:)>(2/3)*Laxon & backwardi(1,2,:)<Laxon));
percentbackdisi=(nbackwarddisi/(nforwarddisi+nbackwarddisi))*100;

%Middle region, defined as Laxon/3<x<2*Laxon/3
nforwardmidi=numel(forwardi(1,2,forwardi(1,2,:)>(1/3)*Laxon & forwardi(1,2,:)<(2/3)*Laxon));
nbackwardmidi=numel(backwardi(1,2,backwardi(1,2,:)>(1/3)*Laxon & backwardi(1,2,:)<(2/3)*Laxon));
percentbackmidi=(nbackwardmidi/(nforwardmidi+nbackwardmidi))*100;

% Overall
nbackwardtoti=numel(backwardi(1,2,backwardi(1,2,:)>0 & backwardi(1,2,:)<Laxon));
nforwardtoti=numel(forwardi(1,2,forwardi(1,2,:)>0 & forwardi(1,2,:)<Laxon));
percenttoti=(nbackwardtoti/(nforwardtoti+nbackwardtoti))*100;

figure
xbar=categorical({'Proximal' 'Middle' 'Distal' 'Overall'});
xbar=reordercats(xbar,{'Proximal' 'Middle' 'Distal' 'Overall'});
ybar=[percentbackproxi percentbackmidi percentbackdisi percenttoti];
bar(xbar,ybar)
ylabel('Percent retrograde MTs')
title('Initial polarity distribution')
    
    
% Program calculation loop
Polarity=zeros(Laxon+1,3); % Initialize polarity array. First column: number of plus-out, 
                            %2nd column, number of minus out; 3rd column, fraction of minus out 
for i=1:ntime
    
    % update length and time for i+1
    A(i+1,1,:)=A(i,1,:)+dt;
    %A(i+1,9,:)=A(i,9,:); % No MTs change length (for this version of program)
    
    % Polarity pattern
      for loc=1:Laxon
          for j=1:NMT
              if (round(A(i,2,j)+1==loc || round(A(i,2,j)+A(i,10,j))+1==loc ...
                || (round(A(i,2,j)+1<loc & round(A(i,2,j)+A(i,10,j))+1>loc)))) % MT overlaps region
                  if(A(i,11,j)==0) % plus-end-out
                  Polarity(loc,1)=Polarity(loc,1)+1; %Number of plus-out MTs at each location
                  elseif(A(i,11,j)==1) %minus-end-out
                  Polarity(loc,2)= Polarity(loc,2)+1; % Number of minus-out MTs at each location
                  end
              end
          end
          if(Polarity(loc,1)+Polarity(loc,2)>0)
          Polarity(loc,3)=Polarity(loc,2)/(Polarity(loc,1)+Polarity(loc,2));
          end
      end
     
     
    
    % Placeholder arrays for attachment numbers
    Nd1=A(i,6,:);
    Nd2=A(i,7,:);
    Nd3=A(i,8,:);
    Nk1=A(i,9,:);
    Nk2=A(i,10,:);
    Nk3=A(i,11,:);
    Nxpar=A(i,12,:);
    Nxanti=A(i,13,:);
    
    Np=Nd1+Nd3+Nk2; % Sum of motors pulling MT with its plus end leading
    Nm=Nd2+Nk1+Nk3; % Sum of motors pulling MT with its minus end leading
    
    
    %Allow motors and crosslinkers to attach where MTs overlap
    for j=1:NMT
        
        % Structure of j loop: (1) Update attachment numbers for time j+1 based on variables
        % at time j (which were determined in the last loop). (2) With
        % updated attachment numbers for time j+1, calculate new variables
        % for j+1, assign them to the appropriate cell.
        
       
        
        % Update position at time i+1 based on position and velocity at time i
        
        A(i+1,2,j)=A(i,2,j)+A(i,3,j)*dt; % Translocate MT based on vel calculated for this time step
        A(i+1,15,j)=A(i,15,j); % Default is for orientation to stay the same. Several possibilities later for this to flip.
        % Apply conditions for if MT reaches a boundary
        if (A(i+1,2,j)<0) % MT has been cleared from axon into cell body
            'MT cleared'
            % Replaced cleared MT with a new short growing MT, random
            % location and orientation
            A(i+1,2,j)=Laxon*rand([1,1,1]); % random location
            A(i+1,3,j)=0.0; % Initial velocity (v)
            A(i+1,4,j)=0.0; % Initial force (Fdp)
            A(i+1,5,j)=0.0; % Initial force (Fdm)
            A(i+1,6,j)=0.0; % Initial Nd1
            A(i+1,7,j)=0.0; % Initial Nd2
            A(i+1,8,j)=0.0; % Initial Nd3
            A(i+1,9,j)=0.0; % Initial Nk1
            A(i+1,10,j)=0.0; % Initial Nk2
            A(i+1,11,j)=0.0; % Initial Nk3
            A(i+1,12,j)=0.0; % Initial Nxpar
            A(i+1,13,j)=0.0; % Initial Nxanti
            A(i+1,14,j)=0.1; % Short initial length, 0.1micron for newly nucleated MT
            
            % Fraction of plus-out MTs at new MT's location
            loc=round(A(i+1,2,j));
            Fmin=Polarity(loc+1,3);
            if(rand<Fmin) % Newly nucleated MT has random orientation
                A(i+1,15,j)=1; % minus-end-out
                'minus out new MT'
            else
                A(i+1,15,j)=0; % plus-end-out
                'plus out new MT'
            end
            
            
            A(i+1,16,j)=1; % dynamic and growing
        elseif(A(i+1,2,j)>Laxon) % MT hits distal end
            A(i+1,2,j)=Laxon; % Can't grow further
            A(i+1,16,j)=0; % Switches to stable
        end
        
        % Update length of dynamic MTs
        if(A(i,16,j)==0)
            A(i+1,14,j)=A(i,14,j);
        elseif(A(i,16,j)==1)
            A(i+1,14,j)=A(i,14,j)+vp*dt;
            if(A(i+1,14,j)>Laxon)
                A(i+1,14,j)=Laxon;
            end
        elseif(A(i,16,j)==2)
            A(i+1,14,j)=A(i,14,j)-vm*dt;
        end
        
        % Apply conditions for if MT length shrinks to zero, nucleate new MT
        if(A(i+1,14,j)<0)
           
            A(i+1,3,j)=0.0; % Initial velocity (v)
            A(i+1,4,j)=0.0; % Initial force (Fdp)
            A(i+1,5,j)=0.0; % Initial force (Fdm)
            A(i+1,6,j)=0.0; % Initial Nd1
            A(i+1,7,j)=0.0; % Initial Nd2
            A(i+1,8,j)=0.0; % Initial Nd3
            A(i+1,9,j)=0.0; % Initial Nk1
            A(i+1,10,j)=0.0; % Initial Nk2
            A(i+1,11,j)=0.0; % Initial Nk3
            A(i+1,12,j)=0.0; % Initial Nxpar
            A(i+1,13,j)=0.0; % Initial Nxanti
            A(i+1,14,j)=0.1; % Short initial length, 0.1micron for newly nucleated MT
            
            % 
            A(i+1,2,j)=Laxon*rand([1,1,1]); % random location
            
               % Fraction of plus-out MTs at new MT's location
            loc=round(A(i+1,2,j));
            Fmin=Polarity(loc+1,3);
            if(rand<Fmin) % Newly nucleated MT has random orientation
                A(i+1,15,j)=1; % minus-end-out
                %'minus out new MT'
            else
                A(i+1,15,j)=0; % plus-end-out
                %'plus out new MT'
            end
            
            A(i+1,16,j)=1; % dynamic and growing 
        end
        
        
        % Give MTs a chance to switch between stable and dynamic
        % 0 = stable, 1 = dynamic & growing, 2= dynamic and shrinking
        if(rand<rdyn*dt & A(i,16,j)==0)
            if(rand<0.5)
                A(i+1,16,j)=1;
            else
                A(i+1,16,j)=2;
            end
        elseif(rand<fpm*dt & A(i,16,j)==1)
            A(i+1,16,j)=2;
            'catastrophe'
        elseif(rand<fmp*dt & A(i,16,j)==2)
            A(i+1,16,j)=1;
            'rescue'
        else
            A(i+1,16,j)=A(i,16,j);
        end
           
        
         % Abbreviated variable names for MT j at time i 
        t=A(i,1,j);
        x=A(i,2,j);
        v=A(i,3,j);
        Fp=A(i,4,j);
        Fm=A(i,5,j);
        L=A(i,14,j);
        P=A(i,15,j);
        
        
        % (1) Update attachment numbers
        %---------------------------------------------------------
        %% Detachment events
        % Calculate detachment rates for motors and non-motile cross-linkers
        
        % Detachment rate for Nd1, Nd3, and Nk2 ('forward pulling' motors)
        
        if(Fp<Np*Fs)
            roffd1=(doff*exp(Fp/(Np(j)*Fdd))+2*v/L)*dt; % Substall detachment rate for Nd1
            roffd3=doff*exp(Fp/(Np(j)*Fdd))*dt; % Substall detachment rate for Nd3
            roffk2=koff*exp(Fp/(Np(j)*Fdk))*dt; % Substall detachment for Nk2 
        else
            roffd1=(doff*(4/(1-exp(-Fp/(Np(j)*1.97))))+2*v/L)*dt; % Superstall detachment rate for Nd1
            roffd3=doff*(4/(1-exp(-Fp/(Np(j)*1.97))))*dt; % Superstall detachment rate for Nd3
            roffk2=koff*(1.54+0.19*(Fp/Np(j)))*dt; % Superstall detachment for Nk2 
        end
        
        % Detachement rate for Nd2, Nk1, Nk3 (backward pulling motors)
        
        if(Fm<Nm*Fs)
            roffd2=doff*exp(Fm/(Nm(j)*Fdd))*dt; % Substall detachment rate for Nd2
            roffk1=(koff*(exp(Fm/(Nm(j)*Fdk)))+2*v/L)*dt; % Substall detachment rate for Nk1 
            roffk3=koff*(exp(Fm/(Nm(j)*Fdk)))*dt; % Substall detachment rate for Nk3 
        else
            roffd2=doff*(4/(1-exp(-Fm/(Nm(j)*1.97))))*dt; % Superstall detachment rate for Nd2
            roffk1=(koff*(1.54+0.19*(Fm/Nm(j)))+2*v/L)*dt; % Superstall detachment rate for Nk1 (update with function)
            roffk3=koff*(1.54+0.19*(Fm/Nm(j)))*dt; % Superstall detachment rate for Nk3 (update with function)
        end
        
        % Calculate detachment rates rates for cross-linkers
        
        roffxpar=xoffpar*dt; % Detachment rate for cross-linkers between parallel MTs
        roffxanti=xoffanti*dt; % Detachment rate for cross-linkers between anti-parallel MTs
        
        % Allow detachment events. Assign to variable name as a placeholder. Assign new values to j+1 after all
        % possible attachment/detachment possibilities have occured for this time step.
        
        % Allow plus-pulling motor detachment events
        
        if(rand<roffd1)
            Nd1(j)=Nd1(j)-1; 
            if(Nd1(j)<0) 
                Nd1(j)=0;
            end
        end
        
        if(rand<roffd2)
            Nd2(j)=Nd2(j)-1;
            if(Nd2(j)<0)
                Nd2(j)=0;
            end
        end
        
        if(rand<roffd3)
            Nd3(j)=Nd3(j)-1;
            if(Nd3(j)<0)
                Nd3(j)=0;
            end
        end
        
        if(rand<roffk1)
            Nk1(j)=Nk1(j)-1;
            if(Nk1(j)<0)
                Nk1(j)=0;
            end
        end
        
        if(rand<roffk2)
            Nk2(j)=Nk2(j)-1;
            if(Nk2(j)<0)
                Nk2(j)=0;
            end
        end
        
        if(rand<roffk3)
            Nk3(j)=Nk3(j)-1;
            if(Nk3(j)<0)
                Nk3(j)=0;
            end
        end
        
        if(rand<roffxpar)
            Nxpar(j)=Nxpar(j)-1;
            if(Nxpar(j)<0)
                Nxpar(j)=0;
            end
        end
        
        if(rand<roffxanti)
            Nxanti(j)=Nxanti(j)-1;
            if(Nxanti(j)<0)
                Nxanti(j)=0;
            end
        end
        
         
       %% Attachment events     
        % Calculate attachment rates
        
        % Find length of overlap between MTs j and k
        for k=j+1:NMT
  
            %Find overlap between MT j and MT k
             if(A(i,2,k)>A(i,2,j)) % If lhs of MT k is to right of lhs of MT j
                 overlap=A(i,2,j)+A(1,14,j)-A(i,2,k); % rhs of MT j - lhs of MT k (length of their overlap)
                 if(overlap<0)
                     overlap=0;
                 end
             else % if lhs of MT j is to right of lhs of MT k
               overlap=A(i,2,k)+A(1,14,k)-A(i,2,j); % rhs of MT k - lhs of MT j (length of overlap between them)
                 if(overlap<0)
                     overlap=0;
                 end
             end
            
             % total number of binding sits on MTs j and k
             Ntotj=round(Dmax*A(i,14,j)); % Number of binding sites on MT j
             Ntotk=round(Dmax*A(i,14,k)); % Number of binding sites on MT k
             
             % Estimated number of overlapping binding sites for j and k
             
             Noverlapj=(overlap/A(i,14,j))*(Ntotj-(A(i,6,j)+A(i,7,j)+A(i,8,j)+A(i,9,j)+A(i,10,j)+A(i,11,j)+A(i,12,j)+A(i,13,j))); % Number of available sites in j overlapping region
             Noverlapk=(overlap/A(i,10,k))*(Ntotk-(A(i,6,k)+A(i,7,k)+A(i,8,k)+A(i,9,k)+A(i,10,k)+A(i,11,k)+A(i,12,k)+A(i,13,k))); % Number of available sites in k overlapping region
             Noverlap=min([Noverlapj,Noverlapk]); % Smaller number of available sites in overlapping region. Use this to estimate number of possible cross-link events that could take place between j and k.
             if(Noverlap<0)
                 Noverlap=0;
             end
             
             
             % Calculate on-rate for each population for MTs j and k
             % Ron depends if they are parallel or anti-parallel
             
             if(A(i,15,j)==A(i,15,k)) % if j and k are parallel
                 ronxpar = Noverlap*xonpar*dt; % On-rate for cross-linkers between parallel MTs
                 ronxanti=0.0;
             else % j and k are anti-parallel
                 ronxanti = Noverlap*xonantipar*dt; % On-rate for cross-linkers between anti-parallel MTs
                 ronxpar=0.0;
             end
             
             rond = Noverlap*don*dt; % On-rate for dynein motors
             ronk = Noverlap*kon*dt; % On-rate for kinesin motors
             
             % Allow motors and cross-linkers to attach to MTs j and k
             
             if(rand<rond) % dynein cross-link forms between j and k
                 if(rand<0.5) % flip a coin to decide if j gets motor domains
                     Nd1(j)=Nd1(j)+1;
                     if(A(i,15,j)==A(i,15,k)) % MTs parallel
                         Nd2(k)=Nd2(k)+1;
                     else
                         Nd3(k)=Nd3(k)+1;
                     end
                 else
                     Nd1(k)=Nd1(k)+1;
                     if(A(i,15,j)==A(i,15,k))
                         Nd2(j)=Nd2(j)+1;
                     else
                         Nd3(j)=Nd3(j)+1;
                     end
                 end
             end
             
             if(rand<ronk) % kinesin cross-link forms between j and k
                 if(rand<0.5) % flip a coin to decide if j gets motor domains
                     Nk1(j)=Nk1(j)+1;
                     if(A(i,15,j)==A(i,15,k)) % MTs parallel
                         Nk2(k)=Nk2(k)+1;
                     else
                         Nk3(k)=Nk3(k)+1;
                     end
                 else
                     Nk1(k)=Nk1(k)+1;
                     if(A(i,15,j)==A(i,15,k))
                         Nk2(j)=Nk2(j)+1;
                     else
                         Nk3(j)=Nk3(j)+1;
                     end
                 end
             end
              
             if (rand < ronxpar) % if random number < parallel-cross-linker rate, then both MTs gain a parallel x-linker
                 Nxpar(j)=Nxpar(j)+1;
                 Nxpar(k)=Nxpar(k)+1;
             end
             
             if (rand < ronxanti) % if random number < anti-parallel-cross-linker rate, then both MTs gain an anti-parallel x-linker
                 Nxanti(j)=Nxanti(j)+1;
                 Nxanti(k)=Nxanti(k)+1;
             end
             
        end % end loop over j - k MT combos, which was used for adding new motors and cross-linkers between two MTs
        
         % (3) Assign updated attachment numbers to i+1
     A(i+1,6,:)=Nd1(:);
     A(i+1,7,:)=Nd2(:);
     A(i+1,8,:)=Nd3(:);
     A(i+1,9,:)=Nk1(:);
     A(i+1,10,:)=Nk2(:);
     A(i+1,11,:)=Nk3(:);
     A(i+1,12,:)=Nxpar(:);
     A(i+1,13,:)=Nxanti(:);
     
      Np=Nd1+Nd3+Nk2; % Sum of motors pulling MT with its plus end leading
      Nm=Nd2+Nk1+Nk3; % Sum of motors pulling MT with its minus end leading
      
      %%
         
        % (2) Allow short MTs to flip orientation (perform operation for
        % i+1 values since they have all been updated at this point. This
        % prevents overwriting orientation update that occured when MTs
        % were cleared from axon or depolymerized to nothing. 
        % ----------------------------------------------------------

             
         if (A(i+1,14,j)<2.0 & Np(j)==0 & Nm(j)==0 & Nxpar(j)==0 & Nxanti(j)==0) % MT short and unattached (length < 2microns and all motor and cross-linkers attachment numbers = 0)
            if(rand<rflip*dt) % if random number < flipping probability, then flipping event occurs
                if (A(i+1,15,j)==0) 
                    A(i+1,15,j)=1; %plus-end-out MT becomes minus-end-out 
                    'flip to minus'
                else
                    A(i+1,15,j)=0; % minus-end-out MT becomes plus-end-out
                    'flip to plus'
                end
           % else
           %     A(i+1,11,j)=A(i,11,j); % else, MT keeps its original orientation
            end
        % else %else if MT *not* short and unattached, it keeps its original orientation
         %    A(i+1,11,j)=A(i,11,j);
         end
    
  %%

     
     % (4) Update v, Fdp for time i+1 based on new i+1 attachment numbers
      if (Np(j) >= Nm(j)) % if MT j has more plus-pulling dynein than minus-pulling dynein (plus team 'winning') (or equal)
       
          %update v
        if(A(i+1,15,j)==0) % if MT j has plus-end-out orientation
            A(i+1,3,j)=vf*(Np(j)-Nm(j))./(Np(j)+(vf/vb)*Nm(j)+...
          (vf/Fs)*(xi+Nxpar(j)*gammapar+Nxanti(j)*gammaanti)); % Calculate new velocity (this equation for velocity is derived from balance of forces)
        elseif(A(i+1,15,j)==1) % if MT j has minus-end-out orientation
            A(i+1,3,j)=-vf*(Np(j)-Nm(j))./(Np(j)+(vf/vb)*Nm(j)+...
            (vf/Fs)*(xi+Nxpar(j)*gammapar+Nxanti(j)*gammaanti)); % Calculate new velocity (opposite direction than plus-end-out MT)
        end
        
        % Update Fdp (force from dynein motors pulling MT j in plus direction)
        Fnum=(Np(j).*(Nm(j)*Fs*(1+(vf/vb))+...
          vf*Np(j).*(xi+Nxpar(j)*gammapar+Nxanti(j)*gammaanti))); % numerator of force expression (equation derived based on force balance)
      
        Fdenom=Np(j)+(vf/vb)*Nm(j)+...
          (vf/Fs)*(xi+Nxpar(j)*gammapar+Nxanti(j)*gammaanti); % denominator of force expression
      
        A(i+1,4,j)= Fnum/Fdenom; % full expression for updated Fdp (force from dynein pulling MT j in plus direction)
        A(i+1,5,j)= A(i+1,4,j)-((Nxpar(j)*gammapar+Nxanti(j)*gammaanti+xi)*A(i+1,3,j)); % Update Fdm (force from dynein pulling MT j in minus direction)
        
      elseif (Nm(j)>Np(j)) % if MT j has more minus-pulling dynein than plus-pulling dynein (minus team 'winning')
          % Update v
          if (A(i+1,15,j)==0) % if MT j has plus-end-out polarity
              A(i+1,3,j)=-vf*(Nm(j)-Np(j))./(Nm(j)+(vf/vb)*Np(j)+...
                  (vf/Fs)*(xi+Nxpar(j)*gammapar+Nxanti(j)*gammaanti)); % new velocity
          elseif(A(i+1,15,j)==1) % if MT j has minus-end-out polarity
              A(i+1,3,j)=vf*(Nm(j)-Np(j))./(Nm(j)+(vf/vb)*Np(j)+...
                  (vf/Fs)*(xi+Nxpar(j)*gammapar+Nxanti(j)*gammaanti)); % new velocity (opposite direction than plus-end-out MT)
          end
          
          % Update Fdp (force from dynein pulling MT j in plus direction)
          
          Fnum=(1+(vf/vb))*Np(j)*Nm(j)*Fs+...
              vf*Nm(j)*(xi+Nxpar(j)*gammapar+Nxanti(j)*gammaanti); % numerator of force expression
          Fdenom=Nm(j)+(vf/vb)*Np(j)+...
          (vf/Fs)*(xi+Nxpar(j)*gammapar+Nxanti(j)*gammaanti); % denominator of force expression
          
          A(i+1,5,j)=Fnum/Fdenom; % full expression for updated Fdp (force from dynein pulling MT j in plus direction)
          A(i+1,4,j)=A(i+1,5,j)+((Nxpar(j)*gammapar+Nxanti(j)*gammaanti+xi)*A(i+1,3,j)); % Update Fdm (force from dynein pulling MT j in minus direction)
         
      end % end if statement for updating force and velocity for MT j at time i+1
    
    end % end loop over different MTs   
end % end loop over time steps             
          
% Graph and save output   
     
% Plot figures of the results

 figure 
 plot([0:Laxon]',Polarity(:,3))
 title('Fraction of minus-end-out')

figure
plot(A(:,1,1),A(:,2,1),A(:,1,2),A(:,2,2),A(:,1,3),A(:,2,3),A(:,1,4),A(:,2,4),A(:,1,5),A(:,2,5),A(:,1,6),A(:,2,6))
xlabel('Time (s)')
ylabel('Position (\mum)')
title('Microtubule position vs. time, sample trajectories')
legend('MT 1','MT 2','MT 3','MT 4','MT 5','MT 6')

figure
plot(A(:,1,1),A(:,6,1),A(:,1,1),A(:,7,1),A(:,1,1),A(:,8,1),A(:,1,1),A(:,9,1),A(:,1,1),A(:,10,1),A(:,1,1),A(:,11,1),A(:,1,1),A(:,12,1),A(:,1,1),A(:,13,1))
xlabel('Time (s)')
ylabel('Attachment number')
legend('N_{d1}','N_{d2}','N_{d3}','N_{k1}','N_{k2}','N_{k3}','Parallel MT Cross-linkers','Anti-parallel MT Cross-linkers')
title('Motors and cross-linkers attached to MT 1')

% Make a video of polarity distribution over time
% Initialize video
myVideo = VideoWriter('MT_t_060221'); %open video file (change file name to avoid overwriting)
myVideo.FrameRate = 10;  %can adjust this, 5 - 10 works well for me
open(myVideo)

nmov=200; %number of time steps to include in movie (ntime must be multiple of nmov)
figure
for j=1:nmov
    tnow=1+(ntime/nmov)*(j);
    %tnow=1+10*j;
    forward=A(tnow,:,A(tnow,15,:)==0);
    backward=A(tnow,:,A(tnow,15,:)==1);
    histogram(forward(1,2,:))
    hold on
    histogram(backward(1,2,:))
    xlim([0 Laxon])
    ylim([0 100])
    ylabel('Total number of MTs')
    xlabel('Distance from cell body, \mum')
    legend('plus-end-out','minus-end-out')
    hold off
    pause(0.1)
    frame = getframe(gcf); %get frame
    writeVideo(myVideo, frame);
      
end
close(myVideo)

% Calculate percent of + out and - out MTs in each region for bar plot

%Proximal region, defined as 0<x<Laxon/3
nforwardprox=numel(forward(1,2,forward(1,2,:)>0 & forward(1,2,:)<Laxon/3)); % # plus-out MTs in this region
nbackwardprox=numel(backward(1,2,backward(1,2,:)>0 & backward(1,2,:)<Laxon/3)); % # of minus-out MTs in this region
percentbackprox=(nbackwardprox/(nforwardprox+nbackwardprox))*100; % percent with minus-out

% Distal region, defined as x>(2/3)*Laxon
nforwarddis=numel(forward(1,2,forward(1,2,:)>(2/3)*Laxon & forward(1,2,:)<Laxon));
nbackwarddis=numel(backward(1,2,backward(1,2,:)>(2/3)*Laxon & backward(1,2,:)<Laxon));
percentbackdis=(nbackwarddis/(nforwarddis+nbackwarddis))*100;

%Middle region, defined as Laxon/3<x<2*Laxon/3
nforwardmid=numel(forward(1,2,forward(1,2,:)>(1/3)*Laxon & forward(1,2,:)<(2/3)*Laxon));
nbackwardmid=numel(backward(1,2,backward(1,2,:)>(1/3)*Laxon & backward(1,2,:)<(2/3)*Laxon));
percentbackmid=(nbackwardmid/(nforwardmid+nbackwardmid))*100;

% Overall
nbackwardtot=numel(backward(1,2,backward(1,2,:)>0 & backward(1,2,:)<Laxon));
nforwardtot=numel(forward(1,2,forward(1,2,:)>0 & forward(1,2,:)<Laxon));
percenttot=(nbackwardtot/(nforwardtot+nbackwardtot))*100;

figure
xbar=categorical({'Proximal' 'Middle' 'Distal' 'Overall'});
xbar=reordercats(xbar,{'Proximal' 'Middle' 'Distal' 'Overall'});
ybar=[percentbackprox percentbackmid percentbackdis percenttot];
bar(xbar,ybar)
ylabel('Percent retrograde MTs')
title('Final polarity distribution')

figure
scatter(A(2,14,A(2,15,:)==0),A(2,3,A(2,15,:)==0))
hold on
scatter(A(2,14,A(2,15,:)==1),A(2,3,A(2,15,:)==1))
legend('plus-end-out','minus-end-out')
xlabel('MT length (\mum)')
ylabel('MT velocity (\mum/min)')
title('Initial velocity vs. length')
hold off


figure
scatter(A(ntime,14,A(ntime,15,:)==0),A(ntime,3,A(ntime,15,:)==0))
hold on
scatter(A(ntime,14,A(ntime,15,:)==1),A(ntime,3,A(ntime,15,:)==1))
legend('plus-end-out','minus-end-out')
xlabel('MT length (\mum)')
ylabel('MT velocity (\mum/min)')
title('Final velocity vs. length')
hold off

figure
scatter(A(2,2,A(2,15,:)==0),A(2,3,A(2,15,:)==0))
hold on
scatter(A(2,2,A(2,15,:)==1),A(2,3,A(2,15,:)==1))
legend('plus-end-out','minus-end-out')
xlabel('Location of MT (\mum)')
ylabel('Velocity of MT (\mum/min)')
title('Initial velocity vs. location')
hold off

figure
scatter(A(ntime,2,A(ntime,15,:)==0),A(ntime,3,A(ntime,15,:)==0))
hold on
scatter(A(ntime,2,A(ntime,15,:)==1),A(ntime,3,A(ntime,15,:)==1))
legend('plus-end-out','minus-end-out')
xlabel('Location of MT (\mum)')
ylabel('Velocity of MT (\mum/sec)')
title('Final velocity vs. location')
hold off

figure
histogram(A(ntime,3,A(ntime,15,:)==0))
hold on
histogram(A(ntime,3,A(ntime,15,:)==1))
legend('plus-end-out','minus-end-out')
ylabel('Count')
xlabel('Velocity of MT (\mum/sec)')
title('Final velocities')
hold off

figure
histogram(A(ntime,6,A(ntime,15,:)==0))
hold on
histogram(A(ntime,7,A(ntime,15,:)==0))
histogram(A(ntime,8,A(ntime,15,:)==0))
title('Histogram of dynein attachment numbers for plus-out MTs')
xlabel('Attachment number')
ylabel('Count')
legend('N_{d1} (sliding filament, +F)','N_{d2} (walking on parallel MT, -F)','N_{d3} (walking on anti-parallel MT, +F')
hold off

figure
histogram(A(ntime,6,A(ntime,15,:)==1))
hold on
histogram(A(ntime,7,A(ntime,15,:)==1))
histogram(A(ntime,8,A(ntime,15,:)==1))
title('Histogram of dynein attachment numbers for minus-out MTs')
xlabel('Attachment number')
ylabel('Count')
legend('N_{d1} (sliding filament, -F)','N_{d2} (walking on parallel MT, +F)','N_{d3} (walking on anti-parallel MT, -F')
% Write out data (writes out LOTS of files, so uncomment when sure you
% want to save the data). Change filename to avoid overwriting

%for b=1:NMT
%temp=A(:,:,b);
%save(sprintf('A_5_101220_%d.txt',b),'temp','-ASCII','-double');
%end

% Write out parameter file (change filename to avoid overwriting)
fid=fopen('A_t_060921_parms.txt','wt');
fprintf(fid,'Parameters, A_t_060921\r\n');
fprintf(fid,'-----------------------\r\n');
fprintf(fid,'Program: Axonal_MTs_dyn_kin_060921.m\r\n');
fprintf(fid,'don\t%f\r\n',don);
fprintf(fid,'kon\t%f\r\n',kon);
fprintf(fid,'xonpar\t%f\r\n',xonpar);
fprintf(fid,'xonantipar\t%f\r\n',xonantipar);
fprintf(fid,'xoffpar\t%f\r\n',xoffpar);
fprintf(fid,'xoffanti\t%f\r\n',xoffanti);
fprintf(fid,'gammapar\t%f\r\n',gammapar);
fprintf(fid,'gammaanti\t%f\r\n',gammaanti);
fprintf(fid,'pflipi\t%f\r\n',pflipi);
fprintf(fid,'rflip\t%f\r\n',rflip);
fclose(fid);


 