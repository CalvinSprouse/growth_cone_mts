delta = 0.0001; % time step
N = 40;  % number of grid points
L = 12.00; % the system extends from x=0 to x=L
deltax = L/(N-1);  x=deltax*(0:N-1);% grid size
h=delta/deltax; hh=delta/(deltax^2); nstep=100/delta; 
%nstep=50;

lambda=3; F=4; gamma=1000; b=15; Vstar=0.3;
%v0=20;mu0=3;F0=5700;T0=500;Vstar=6/v0;moff=0.01;beta=0;

v=(N-(1:N))/(N-1);v_new=(N-(1:N))/(N-1); % vel
Fmyo=ones(size(1:N));T=(N-(1:N));
m=ones(1,N); 
m_new = ones(1,N); % myosin distribution


%MAIN LOOP

for j=1:nstep

        v_new(2:N-1) = v(2:N-1) + hh*(v(3:N) + v(1:N-2)...
        -2*v(2:N-1))+(delta/(lambda))*(-T(2:N-1)+Fmyo(2:N-1));
    
        v_new(N) = v(N-1);
        v = v_new;
        
            m_new(2:(N))=m(2:(N))+delta*(1-m(2:(N)))-...
            h*gamma*(v(2:(N)).*m(2:(N))-v(1:(N-1)).*m(1:(N-1)));
            m_new(1)=0;
            m=m_new;
       
    % update myosin stress
         delm(2:N)=(m(2:(N))-m(1:(N-1)))/deltax;
         delm(1)=(m(2)-m(1))/deltax;
         Fmyo=F*m;
        
         Nadh=(1+b*Fmyo).*exp(-v/Vstar);
        T=Nadh.*v;       
        
end


for i=1:40
    if(T(i)==max(T))
        maxi=i;
    end
end
xmax=deltax*(maxi-1)
vmax=v(maxi)

%OUTPUT
x=deltax*(0:(N-1));
Tv=v.*exp(-v/Vstar);
Tf=1+b*Fmyo;

figure
plot(x,v,'m',x,Fmyo,'k--',...
    x,T,'r','LineWidth',2)
legend('Actin flow speed','Myosin density','Traction')

A=horzcat(x',v',T',m',Tv',Tf',Nadh');
Parms=vertcat(lambda,F,gamma,b,Vstar,xmax,vmax);
%Parms=vertcat(v0,mu0,F0,T0,Vstar,moff,beta,xmax,vmax);

dlmwrite('Traction_15_052113.txt',A,'\t')
dlmwrite('Traction_15_parms_052113.txt',Parms)