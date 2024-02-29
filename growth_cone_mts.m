% growth_cone_mts.m

%{
Script description:
This script is a modification to 'Traction_model_052813_3.m' to include
microtubule (mt) action in the growth cone (gc). See 'Craig et al. 2015'
for more information on original model features.

Model description:


Model questions for Dr. Craig:

%}

% configure simulation parameters
% define the system extent or length [micrometer]
L = 12;
% define the number of grid points [ ]
N = 40;
% define the time step [s]
Deltat = 0.0001;

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


% MAIN LOOP
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


% some yet unkown plotting preliminary operation
for i=1:40
    if(T(i)==max(T))
        maxi=i;
    end
end
xmax = Deltax*(maxi-1)
vmax = v(maxi)

% OUTPUT
x = Deltax*(0:(N-1));
Tv = v.*exp(-v/Vstar);
Tf = 1+b*Fmyo;

figure
plot(x,v,'m',x,Fmyo,'k--',x,T,'r','LineWidth',2)
legend('Actin flow speed','Myosin density','Traction')

% save simulation output to .mat file
A = horzcat(x',v',T',M',Tv',Tf',Nadh');
Parms = vertcat(lambda,F,gamma,b,Vstar,xmax,vmax);
% Parms = vertcat(v0,mu0,F0,T0,Vstar,moff,beta,xmax,vmax);
save("growth_cone_mts", "-v7", "-nocompression");
