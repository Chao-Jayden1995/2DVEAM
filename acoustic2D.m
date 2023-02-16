%%%% 2D viscoacoustic wave modeling%%%%%%%%%%%%%%

clear; close all; clc;
isnap=500;Q=100;N=3;

%Parsing parameters
nt       =  2500;
nx       =  1000;nz       =  1000;
dx       =  2;  dz       =  2;

c     =  2000*ones(nz,nx);
rho   =  1000;
K=rho*c.^2;

t_total = 1;                       % [sec] recording duration
dt=4e-04;
nt = round(t_total/dt);             % number of time steps
t = [0:nt]*dt;

f0       =  30;
t0       = 0.1;


jsrc = 500;                 % source location along OZ
isrc = 500;                 % source location along OX
%% ABSORBING BOUNDARY (GSRM)
abs_thick = 120;         % thicknes of the layer
yita(nz,nx)=0;
yita_max=6.3*f0;
Amax=log(1.0/0.01)/2;
Bmax=2;
Dismax=abs_thick*dz;


for iz = 1:abs_thick
   dis=(abs_thick-iz+1)/abs_thick*dz;
   yita(iz,:)=yita_max*dis^Bmax;
   c(iz,:)=c(iz,:)*exp(-Amax*dis^Bmax);
end
for ix = 1:abs_thick
   dis=(abs_thick-ix+1)/abs_thick*dx;
   yita(:,ix)=yita_max*dis^Bmax;
   c(:,ix)=c(:,ix)*exp(-Amax*dis^Bmax);
end
for iz = nz-abs_thick:nz
   dis=(abs_thick-nz+iz)/abs_thick*dz;
   yita(iz,:)=yita_max*dis^Bmax;
   c(iz,:)=c(iz,:)*exp(-Amax*dis^Bmax);
end
for ix = nx-abs_thick:nx
   dis=(abs_thick-nx+ix)/abs_thick*dx;
   yita(:,ix)=yita_max*dis^Bmax;
   c(:,ix)=c(:,ix)*exp(-Amax*dis^Bmax);
end

for iz = 1:abs_thick
for ix = 1:abs_thick    
   dis=sqrt((abs_thick-ix+1)^2+(abs_thick-iz+1)^2)/abs_thick*dx;
   yita(iz,ix)=yita_max*dis^Bmax;
   c(iz,ix)=c(iz,ix)*exp(-Amax*dis^Bmax);
end
end

for iz = 1:abs_thick
for ix = nx-abs_thick:nx   
   dis=sqrt((abs_thick-nx+ix)^2+(abs_thick-iz+1)^2)/abs_thick*dx;
   yita(iz,ix)=yita_max*dis^Bmax;
   c(iz,ix)=c(iz,ix)*exp(-Amax*dis^Bmax);
end
end

for iz =  nz-abs_thick:nz
for ix = 1:abs_thick    
   dis=sqrt((abs_thick-ix+1)^2+(abs_thick-nz+iz)^2)/abs_thick*dx;
   yita(iz,ix)=yita_max*dis^Bmax;
   c(iz,ix)=c(iz,ix)*exp(-Amax*dis^Bmax);
end
end

for iz = nz-abs_thick:nz
for ix = nx-abs_thick:nx   
   dis=sqrt((abs_thick-nx+ix)^2+(abs_thick-nz+iz)^2)/abs_thick*dx;
   yita(iz,ix)=yita_max*dis^Bmax;
   c(iz,ix)=c(iz,ix)*exp(-Amax*dis^Bmax);
end
end

% %%%%%% Qkappa=100


    a=constantQ_GSLS(Q,f0,N);
    
tao_e1=a(1,1);
tao_e2=a(1,2);
tao_e3=a(1,3);

tao_s1=a(2,1);
tao_s2=a(2,2);
tao_s3=a(2,3);



a=pi*(t-t0)*f0;
source=(1-2*a.^2).*exp(-a.^2);


    k = 1;    
    u    = zeros(nz,nx);
    uold = zeros(nz,nx);
    dux  = u; duz = u;
    w1new= zeros(nz,nx);
    w2new= zeros(nz,nx);
    w3new= zeros(nz,nx);
    w1old= zeros(nz,nx);    
    w2old= zeros(nz,nx);   
    w3old= zeros(nz,nx); 
    %####################################################################
    %####                      Begin time loop                       ####
    %####################################################################
    temp=tao_e1/tao_s1+tao_e2/tao_s2+tao_e3/tao_s3;
    cons1=(1-tao_e1/tao_s1)/(tao_s1)*c(iz,ix)^2/temp;
    cons2=(1-tao_e2/tao_s2)/(tao_s2)*c(iz,ix)^2/temp;
    cons3=(1-tao_e3/tao_s3)/(tao_s3)*c(iz,ix)^2/temp;
    e1=exp(-dt/tao_s1);
    e2=exp(-dt/tao_s2);
    e3=exp(-dt/tao_s3);
    %% PLRC
coe12=tao_s1^2*(1.0-e1-e1*dt/tao_s1)/dt;coe11=-tao_s1*(e1-1.0);coe11=coe11-coe12;   
coe22=tao_s2^2*(1.0-e2-e2*dt/tao_s2)/dt;coe21=-tao_s2*(e2-1.0);coe21=coe21-coe22;   
coe32=tao_s3^2*(1.0-e3-e3*dt/tao_s3)/dt;coe31=-tao_s3*(e3-1.0);coe31=coe31-coe32;   
    %% PCRC
% coe11=-tao_s1*(e1-1.0);coe12=0;  
% coe21=-tao_s2*(e2-1.0);coe22=0;
% coe31=-tao_s3*(e3-1.0);coe32=0;
    
    
    
    a1=9.0/5040.0;
    a2=128.0/5040.0;
    a3=1008.0/5040.0;
    a4=8064.0/5040.0;
    a5=-14350.0/5040.0;
    
    p =  zeros(nz,nx);
    
    
    for it=2:nt
        unew = zeros(nz,nx);
        unew(isrc,jsrc) = unew(isrc,jsrc) + dt^2 * source(it);
        
 
         
        for ix=5:nx-4     
            for iz=5:nz-4
                
               uxx=(-a1*u(iz,ix-4)+a2*u(iz,ix-3)-a3*u(iz,ix-2)+a4*u(iz,ix-1)+a5*u(iz,ix)+a4*u(iz,ix+1)-a3*u(iz,ix+2)+a2*u(iz,ix+3)-a1*u(iz,ix+4))/dx^2;
               uzz=(-a1*u(iz-4,ix)+a2*u(iz-3,ix)-a3*u(iz-2,ix)+a4*u(iz-1,ix)+a5*u(iz,ix)+a4*u(iz+1,ix)-a3*u(iz+2,ix)+a2*u(iz+3,ix)-a1*u(iz+4,ix))/dz^2;
     
               
               theta(iz,ix)=uxx+uzz;
%% %%%%%% persure-version acoustic wave equation               
                unew(iz,ix) = unew(iz,ix) + (c(iz,ix)^2) * (dt^2) *(theta(iz,ix)) + 2 * u(iz,ix) - uold(iz,ix)-...
                    yita(iz,ix)*dt*u(iz,ix) +yita(iz,ix)*dt*uold(iz,ix); 

%% %%%%%% persure-version memory variable viscoacoustic wave equation
%              w1new(iz,ix)=(w1old(iz,ix)+dt*(-0.5*w1old(iz,ix)/tao_s1-uxx -uzz))/(1+0.5*dt/tao_s1);
%              w2new(iz,ix)=(w2old(iz,ix)+dt*(-0.5*w2old(iz,ix)/tao_s2-uxx -uzz))/(1+0.5*dt/tao_s2);
%              w3new(iz,ix)=(w3old(iz,ix)+dt*(-0.5*w3old(iz,ix)/tao_s3-uxx -uzz))/(1+0.5*dt/tao_s3);
%              unew(iz,ix) = unew(iz,ix) + (c(iz,ix)^2) * (dt^2) *(uxx +uzz) ...
%               -(dt^2) *(cons1*0.5*(w1new(iz,ix)+w1old(iz,ix))+cons2*0.5*(w2new(iz,ix)+w2old(iz,ix))+cons3*0.5*(w3new(iz,ix)+w3old(iz,ix)))+ 2 * u(iz,ix) - uold(iz,ix)-...
%                     yita(iz,ix)*dt*u(iz,ix) +yita(iz,ix)*dt*uold(iz,ix); ; 

%% %%%%%% persure-version iteration-variable viscoacoustic wave equation
%              uxx_old=(-a1*uold(iz,ix-4)+a2*uold(iz,ix-3)-a3*uold(iz,ix-2)+a4*uold(iz,ix-1)+a5*uold(iz,ix)+a4*uold(iz,ix+1)-a3*uold(iz,ix+2)+a2*uold(iz,ix+3)-a1*uold(iz,ix+4))/dx^2;
%              uzz_old=(-a1*uold(iz-4,ix)+a2*uold(iz-3,ix)-a3*uold(iz-2,ix)+a4*uold(iz-1,ix)+a5*uold(iz,ix)+a4*uold(iz+1,ix)-a3*uold(iz+2,ix)+a2*uold(iz+3,ix)-a1*uold(iz+4,ix))/dz^2;            
%              w1new(iz,ix)=-coe11*(uxx +uzz) - coe12*(uxx_old +uzz_old)+ e1*w1old(iz,ix);
%              w2new(iz,ix)=-coe21*(uxx +uzz) - coe22*(uxx_old +uzz_old)+ e2*w2old(iz,ix);
%              w3new(iz,ix)=-coe31*(uxx +uzz) - coe32*(uxx_old +uzz_old)+ e3*w3old(iz,ix);
%              unew(iz,ix) = unew(iz,ix) + (c(iz,ix)^2) * (dt^2) *(uxx +uzz) - ...
%               (dt^2) *(cons1*w1new(iz,ix)+cons2*w2new(iz,ix)+cons3*w3new(iz,ix))+ 2 * u(iz,ix) - uold(iz,ix)-...
%                     yita(iz,ix)*dt*u(iz,ix) +yita(iz,ix)*dt*uold(iz,ix); ; 
            end
        end


        uold = u  ;
        u    = unew;
        w1old=w1new;
        w2old=w2new;
        w3old=w3new;

        
        
        shot_g(it,1) = it*dt;
        shot_g(it,2) = u(700,700);
       
        
         if (mod(it,isnap)==0)
            U(:,:,k) = u;
            k = k + 1;
        fprintf('Time step: %d \t %.4f s\n',it, single(t(it)));
        imagesc(u); colorbar; colormap jet; xlabel('m'); ylabel('m');
        title(['Step = ',num2str(it),'/',num2str(nt),', Time: ',sprintf('%.4f',t(it)),' sec']);
        axis equal tight; drawnow;    
        end
    end    


