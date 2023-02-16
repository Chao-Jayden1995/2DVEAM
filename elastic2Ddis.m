% --------------------------------------------------------------
clear;clc;
close all;
% Output periodicity in time steps
IT_DISPLAY = 1000;
f0 = 20.0;                          % dominant frequency of the wavelet
t0 = 0.1;                     % excitation time
N=3;

%% SOURCE

factor = 1;                      % amplitude coefficient
angle_force = 180.0;                 % spatial orientation

jsrc = 460;                 % source location along OZ
isrc = 1347;                 % source location along OX
%% MODEL
% Model dimensions, [m]
load('SEG_Hess_Attenuation.mat')
load('SEG_Hess_relaxtion.mat')

[nz,nx]=size(c11);
dx = 2;
dz = 2;
nz=nz-4;
nx=nx-4;
%% TIME STEPPING
t_total = 3.0;                       % [sec] recording duration
dt = 1.5e-4;
nt = round(t_total/dt);             % number of time steps
t = [0:nt]*dt;

%% tao_s; tao_e
rho=1500* ones(nz+4, nx+4);  
temp=(sum(tao_e11./tao_s));a(1,:,:)=c11;temp=a./temp;
for num=1:N
co11(num,:,:)=(1-tao_e11(num,:,:)./tao_s(num,:,:))./(tao_s(num,:,:)).*temp;
end
temp=(sum(tao_e13./tao_s));a(1,:,:)=c13;temp=a./temp;
for num=1:N
co13(num,:,:)=(1-tao_e13(num,:,:)./tao_s(num,:,:))./(tao_s(num,:,:)).*temp;
end
temp=(sum(tao_e33./tao_s));a(1,:,:)=c33;temp=a./temp;
for num=1:N
co33(num,:,:)=(1-tao_e33(num,:,:)./tao_s(num,:,:))./(tao_s(num,:,:)).*temp;
end
temp=(sum(tao_e55./tao_s));a(1,:,:)=c55;temp=a./temp;
for num=1:N
co55(num,:,:)=(1-tao_e55(num,:,:)./tao_s(num,:,:))./(tao_s(num,:,:)).*temp;
end

    %% RCTSE
coe1=tao_s.^2.*(e-1.0+dt./tao_s)/dt;coe2=-tao_s.*(e-1.0);coe2=coe2-coe1;   
    %% PLRC
% coe2=tao_s.^2.*(1.0-e-e.*dt./tao_s)./dt;coe1=-tao_s.*(e-1.0);coe1=coe1-coe2;   
    %% PCRC
% coe1=-tao_s.*(e-1.0);coe2=0.*e./e;  
    %% TRC
% coe1=-0.5.*tao_s.*(e-1.0);coe2=-0.5.*tao_s.*(e-1.0);
    %% RCB
% coe1=0.5.*dt.*(1-dt./tao_s);coe2=0.5.*dt.*e./e;
    %% RCTR
% coe1=0.5*dt.*e./e;coe2=0.5.*dt.*e;



%% ABSORBING BOUNDARY (ABS)
abs_thick =420;         % thicknes of the layer
abs_rate = 0.3/abs_thick;      % decay rate

lmargin = [abs_thick abs_thick];
rmargin = [abs_thick abs_thick];
weights = ones(nz+4,nx+4);
for iz = 1:nz+4
    for ix = 1:nx+4
        i = 0;
        j = 0;
        k = 0;
        if (ix < lmargin(1) + 1)
            i = lmargin(1) + 1 - ix;
        end
        if (iz < lmargin(2) + 1)
            k = lmargin(2) + 1 - iz;
        end
        if (nx - rmargin(1) < ix)
            i = ix - nx + rmargin(1);
        end
        if (nz - rmargin(2) < iz)
            k = iz - nz + rmargin(2);
        end
        if (i == 0 && j == 0 && k == 0)
            continue
        end
        rr = abs_rate * abs_rate * double(i*i + j*j + k*k );
        weights(iz,ix) = exp(-rr);
    end
end



tau=pi*f0*(t-t0);
dt2rho_src = dt/rho(jsrc, isrc);
% source_term = factor * exp(-a*(t-t0).^2);                             % Gaussian
% source_term = factor*(t-t0).*exp(-tau.^2);              % First derivative of a Gaussian:
source_term = factor * (1-2*tau.^2).*exp(-tau.^2);        % Ricker source time function (second derivative of a Gaussian):

force_x = sin(angle_force * pi / 180) * source_term * dt2rho_src / (dx * dz);
force_z = cos(angle_force * pi / 180) * source_term * dt2rho_src / (dx * dz);
min_wavelengh = 0.5*min(min(sqrt(c55./rho)))/f0;     % shortest wavelength bounded by velocity in the air
CFL = max(max(sqrt(c11./rho)))*dt * sqrt(1.0/dx^2 + 1.0/dz^2);
%% SUMMARY
fprintf('#################################################\n');
fprintf('2D elastic FDTD wave propagation in isotripic \nmedium in displacement formulation with \nCerjan(1985) boundary conditions\n');
fprintf('#################################################\n');
fprintf('Model:\n\t%d x %d\tgrid nz x nx\n\t%.1e x %.1e\t[m] dz x dx\n',nz, nx, dz,dx);
fprintf('\t%.1e x %.1e\t[m] model size\n',nx*dx, nz*dz);
% fprintf('\t%.1e...%.1e\t[m/s] vp\n', min(vp(:)), max(vp(:)));
% fprintf('\t%.1e...%.1e\t[m/s] vs\n', min(vs(:)), max(vs(:)));
fprintf('\t%.1e...%.1e\t[kg/m3] rho\n', min(rho(:)), max(rho(:)));
fprintf('Time:\n\t%.1e\t[sec] total\n\t%.1e\tdt\n\t%d\ttime steps\n',t_total,dt,nt);
fprintf('Source:\n\t%.1e\t[Hz] dominant frequency\n\t%.1f\t[sec] index time\n',f0,t0);
fprintf('Other:\n\t%.1f\tCFL number\n', CFL);
fprintf('\t%.2f\t[m] shortest wavelength\n\t%d, %d\t points-per-wavelength OX, OZ\n', min_wavelengh, floor(min_wavelengh/dx), floor(min_wavelengh/dz));
fprintf('#################################################\n');

%% ALLOCATE MEMORY FOR WAVEFIELD
ux3 = zeros(nz+4,nx+4);            % Wavefields at t
uz3 = zeros(nz+4,nx+4);
ux2 = zeros(nz+4,nx+4);            % Wavefields at t-1
uz2 = zeros(nz+4,nx+4);
ux1 = zeros(nz+4,nx+4);            % Wavefields at t-2
uz1 = zeros(nz+4,nx+4);

    dux_dxxo=zeros(nz,nx);
    dux_dxzo=zeros(nz,nx);
    dux_dzzo=zeros(nz,nx);
    duz_dxxo=zeros(nz,nx);
    duz_dxzo=zeros(nz,nx);
    duz_dzzo=zeros(nz,nx);
    dux_dxo=zeros(nz,nx);
    dux_dzo=zeros(nz,nx);
    duz_dxo=zeros(nz,nx);
    duz_dzo=zeros(nz,nx);
    
% Coefficients for derivatives
co_dxx = 1/(12*dx^2);
co_dzz = 1/(12*dz^2);
co_dx = 1/(12*dx);
co_dz = 1/(12*dz);
co_dxz = 1/(12.0 * 12.0 * dx * dz);
co_dzx = 1/(12.0 * 12.0 * dx * dz);
dt2rho=(dt^2)./rho;


c11_dx=zeros(nz,nx);
c13_dx=zeros(nz,nx);
c33_dx=zeros(nz,nx);
c55_dx=zeros(nz,nx);
c11_dz=zeros(nz,nx);
c13_dz=zeros(nz,nx);
c33_dz=zeros(nz,nx);
c55_dz=zeros(nz,nx);
co11_dx=zeros(N,nz,nx);
co13_dx=zeros(N,nz,nx);
co33_dx=zeros(N,nz,nx);
co55_dx=zeros(N,nz,nx);
co11_dz=zeros(N,nz,nx);
co13_dz=zeros(N,nz,nx);
co33_dz=zeros(N,nz,nx);
co55_dz=zeros(N,nz,nx);
tao_s_dx=zeros(N,nz,nx);
tao_s_dz=zeros(N,nz,nx);
ra=zeros(N,nz,nx);
rb=zeros(N,nz,nx);
ra_old=zeros(N,nz,nx);
rb_old=zeros(N,nz,nx);



i_snap=1;
%% Loop over TIME
    % c13
    c13_dx=co_dx * (c13(3:end-2,1:end-4) - 8*c13(3:end-2,2:end-3)+...
        8*c13(3:end-2,4:end-1) - c13(3:end-2,5:end));
    c13_dz=co_dz * (c13(1:end-4,3:end-2) - 8*c13(2:end-3,3:end-2)+...
        8*c13(4:end-1,3:end-2) - c13(5:end,3:end-2));
    % c55
    c55_dx=co_dx * (c55(3:end-2,1:end-4) - 8*c55(3:end-2,2:end-3)+...
        8*c55(3:end-2,4:end-1) - c55(3:end-2,5:end));
    c55_dz=co_dz * (c55(1:end-4,3:end-2) - 8*c55(2:end-3,3:end-2)+...
        8*c55(4:end-1,3:end-2) - c55(5:end,3:end-2));
    % c11
    c11_dx=co_dx * (c11(3:end-2,1:end-4) - 8*c11(3:end-2,2:end-3)+...
        8*c11(3:end-2,4:end-1) - c11(3:end-2,5:end));
    c11_dz=co_dz * (c11(1:end-4,3:end-2) - 8*c11(2:end-3,3:end-2)+...
        8*c11(4:end-1,3:end-2) - c11(5:end,3:end-2));
    % c33
    c33_dx=co_dx * (c33(3:end-2,1:end-4) - 8*c33(3:end-2,2:end-3)+...
        8*c33(3:end-2,4:end-1) - c33(3:end-2,5:end));
    c33_dz=co_dz * (c33(1:end-4,3:end-2) - 8*c33(2:end-3,3:end-2)+...
        8*c33(4:end-1,3:end-2) - c33(5:end,3:end-2));

    for num=1:N
   % co11
    co11_dx(num,:,:)=co_dx * (co11(num,3:end-2,1:end-4) - 8*co11(num,3:end-2,2:end-3)+...
        8*co11(num,3:end-2,4:end-1) - co11(num,3:end-2,5:end));
    co11_dz(num,:,:)=co_dz * (co11(num,1:end-4,3:end-2) - 8*co11(num,2:end-3,3:end-2)+...
        8*co11(num,4:end-1,3:end-2) - co11(num,5:end,3:end-2));
    
    % co13
    co13_dx(num,:,:)=co_dx * (co13(num,3:end-2,1:end-4) - 8*co13(num,3:end-2,2:end-3)+...
        8*co13(num,3:end-2,4:end-1) - co13(num,3:end-2,5:end));
    co13_dz(num,:,:)=co_dz * (co13(num,1:end-4,3:end-2) - 8*co13(num,2:end-3,3:end-2)+...
        8*co13(num,4:end-1,3:end-2) - co13(num,5:end,3:end-2));       
    % co33
    co33_dx(num,:,:)=co_dx * (co33(num,3:end-2,1:end-4) - 8*co33(num,3:end-2,2:end-3)+...
        8*co33(num,3:end-2,4:end-1) - co33(num,3:end-2,5:end));
    co33_dz(num,:,:)=co_dz * (co33(num,1:end-4,3:end-2) - 8*co33(num,2:end-3,3:end-2)+...
        8*co33(num,4:end-1,3:end-2) - co33(num,5:end,3:end-2));       
    % co55
    co55_dx(num,:,:)=co_dx * (co55(num,3:end-2,1:end-4) - 8*co55(num,3:end-2,2:end-3)+...
        8*co55(num,3:end-2,4:end-1) - co55(num,3:end-2,5:end));
    co55_dz(num,:,:)=co_dz * (co55(num,1:end-4,3:end-2) - 8*co55(num,2:end-3,3:end-2)+...
        8*co55(num,4:end-1,3:end-2) - co55(num,5:end,3:end-2));        
    % tao_s
    tao_s_dx(num,:,:)=co_dx * (tao_s(num,3:end-2,1:end-4) - 8*tao_s(num,3:end-2,2:end-3)+...
        8*tao_s(num,3:end-2,4:end-1) - tao_s(num,3:end-2,5:end));
    tao_s_dz(num,:,:)=co_dz * (tao_s(num,1:end-4,3:end-2) - 8*tao_s(num,2:end-3,3:end-2)+...
        8*tao_s(num,4:end-1,3:end-2) - tao_s(num,5:end,3:end-2));
    end
    
    tao_s_dx=tao_s_dx./(tao_s(:,3:end-2,3:end-2).^2);
    tao_s_dz=tao_s_dz./(tao_s(:,3:end-2,3:end-2).^2);
tic;
for it = 1:nt
    ux3 = zeros(size(ux2));
    uz3 = zeros(size(uz2));
    % Second-order derivatives
    % Ux
    dux_dxx = co_dxx * (-ux2(3:end-2,1:end-4) + 16*ux2(3:end-2,2:end-3) - ...
        30*ux2(3:end-2,3:end-2) + 16*ux2(3:end-2,4:end-1)-ux2(3:end-2,5:end));
    dux_dzz = co_dzz * (-ux2(1:end-4,3:end-2) + 16*ux2(2:end-3,3:end-2) - ...
        30*ux2(3:end-2,3:end-2) + 16*ux2(4:end-1,3:end-2)-ux2(5:end,3:end-2));
    dux_dxz = co_dxz * ((ux2(1:end-4,1:end-4) - 8*ux2(1:end-4,2:end-3) + 8*ux2(1:end-4,4:end-1) - ux2(1:end-4,5:end))-...
        8*(ux2(2:end-3,1:end-4) - 8*ux2(2:end-3,2:end-3) + 8*ux2(2:end-3,4:end-1) - ux2(2:end-3,5:end))+...
        8*(ux2(4:end-1,1:end-4) - 8*ux2(4:end-1,2:end-3) + 8*ux2(4:end-1,4:end-1) - ux2(4:end-1,5:end))-...
        (ux2(5:end,1:end-4) - 8*ux2(5:end,2:end-3) + 8*ux2(5:end,4:end-1) - ux2(5:end,5:end)));
    dux_dx=co_dx * (ux2(3:end-2,1:end-4) - 8*ux2(3:end-2,2:end-3)+...
        8*ux2(3:end-2,4:end-1) - ux2(3:end-2,5:end));
    dux_dz=co_dz * (ux2(1:end-4,3:end-2) - 8*ux2(2:end-3,3:end-2)+...
        8*ux2(4:end-1,3:end-2) - ux2(5:end,3:end-2));
    % Uz
    duz_dxx =  co_dxx * (-uz2(3:end-2,1:end-4) + 16*uz2(3:end-2,2:end-3) - ...
        30*uz2(3:end-2,3:end-2) + 16*uz2(3:end-2,4:end-1)-uz2(3:end-2,5:end));
    duz_dzz = co_dzz * (-uz2(1:end-4,3:end-2) + 16*uz2(2:end-3,3:end-2) - ...
        30*uz2(3:end-2,3:end-2) + 16*uz2(4:end-1,3:end-2)-uz2(5:end,3:end-2));
    duz_dxz = co_dxz * ((uz2(1:end-4,1:end-4) - 8*uz2(1:end-4,2:end-3) + 8*uz2(1:end-4,4:end-1) - uz2(1:end-4,5:end))-...
        8*(uz2(2:end-3,1:end-4) - 8*uz2(2:end-3,2:end-3) + 8*uz2(2:end-3,4:end-1) - uz2(2:end-3,5:end))+...
        8*(uz2(4:end-1,1:end-4) - 8*uz2(4:end-1,2:end-3) + 8*uz2(4:end-1,4:end-1) - uz2(4:end-1,5:end))-...
        (uz2(5:end,1:end-4) - 8*uz2(5:end,2:end-3) + 8*uz2(5:end,4:end-1) - uz2(5:end,5:end)));
    duz_dx=co_dx * (uz2(3:end-2,1:end-4) - 8*uz2(3:end-2,2:end-3)+...
        8*uz2(3:end-2,4:end-1) - uz2(3:end-2,5:end));
    duz_dz=co_dz * (uz2(1:end-4,3:end-2) - 8*uz2(2:end-3,3:end-2)+...
        8*uz2(4:end-1,3:end-2) - uz2(5:end,3:end-2));

    

    
%      elastic
%      ux3(3:end-2,3:end-2) = 2.0*ux2(3:end-2,3:end-2) - ux1(3:end-2,3:end-2) ...
%         + (c11(3:end-2,3:end-2) .* dux_dxx + c13(3:end-2,3:end-2) .* duz_dxz +...
%         c55(3:end-2,3:end-2) .* (dux_dzz + duz_dxz) + c11_dx .* dux_dx +...
%         c13_dx .* duz_dz + c55_dz .* (dux_dz + duz_dx)).*dt2rho(3:end-2,3:end-2);
%     
%      uz3(3:end-2,3:end-2) = 2.0*uz2(3:end-2,3:end-2) - uz1(3:end-2,3:end-2) ...
%         + (c55(3:end-2,3:end-2) .* (dux_dxz + duz_dxx) + c13(3:end-2,3:end-2) .* dux_dxz +...
%         c33(3:end-2,3:end-2) .* duz_dzz + c55_dx .* (dux_dz + duz_dx) +...
%         c33_dz .* duz_dz + c13_dz .* dux_dx).*dt2rho(3:end-2,3:end-2);
    %% viscoelastic memory   
   for num=1:N
    ra(num,:,:)=(dt.*(co11(num,3:end-2,3:end-2).*reshape(dux_dxx,[1,nz,nx]) + co13(num,3:end-2,3:end-2).*reshape(duz_dxz,[1,nz,nx]) +...
        co55(num,3:end-2,3:end-2).* reshape((dux_dzz + duz_dxz),[1,nz,nx]) + co11_dx(num,:,:).* reshape(dux_dx,[1,nz,nx])+...
        co13_dx(num,:,:) .* reshape(duz_dz,[1,nz,nx]) + co55_dz(num,:,:) .* reshape((dux_dz + duz_dx),[1,nz,nx]))+(1-0.5*dt./tao_s(num,3:end-2,3:end-2)).*ra_old(num,:,:))./(1+0.5*dt./tao_s(num,3:end-2,3:end-2)); 
   
    rb(num,:,:)=(dt.*(co55(num,3:end-2,3:end-2).*reshape((dux_dxz + duz_dxx),[1,nz,nx]) + co13(num,3:end-2,3:end-2).*reshape(dux_dxz,[1,nz,nx]) +...
        co33(num,3:end-2,3:end-2).*reshape(duz_dzz,[1,nz,nx]) + co55_dx(num,:,:).* reshape((dux_dz + duz_dx),[1,nz,nx])+...
        co33_dz(num,:,:) .* reshape(duz_dz,[1,nz,nx]) + co13_dz(num,:,:) .* reshape(dux_dx,[1,nz,nx]))+(1-0.5*dt./tao_s(num,3:end-2,3:end-2)).*rb_old(num,:,:))./(1+0.5*dt./tao_s(num,3:end-2,3:end-2)); 
  
    end
    
     ux3(3:end-2,3:end-2) = 2.0*ux2(3:end-2,3:end-2) - ux1(3:end-2,3:end-2) + (c11(3:end-2,3:end-2) .* dux_dxx + c13(3:end-2,3:end-2) .* duz_dxz +...
        c55(3:end-2,3:end-2) .* (dux_dzz + duz_dxz) + c11_dx .* dux_dx + c13_dx .* duz_dz + c55_dz .* (dux_dz + duz_dx)+...
        reshape(0.5*sum(ra+ra_old),[nz,nx]) ).*dt2rho(3:end-2,3:end-2);
    
     uz3(3:end-2,3:end-2) = 2.0*uz2(3:end-2,3:end-2) - uz1(3:end-2,3:end-2)+(c55(3:end-2,3:end-2) .* (dux_dxz + duz_dxx) + c13(3:end-2,3:end-2) .* dux_dxz +...
        c33(3:end-2,3:end-2) .* duz_dzz + c55_dx .* (dux_dz + duz_dx) +  c33_dz .* duz_dz + c13_dz .* dux_dx+...
        reshape(0.5*sum(rb+rb_old),[nz,nx]) ).*dt2rho(3:end-2,3:end-2);
    %% viscoelastic iteration    
%  for num=1:N   
%     ra(num,:,:)=coe1(num,3:end-2,3:end-2).*(co11(num,3:end-2,3:end-2).*reshape(dux_dxx,[1,nz,nx]) + co13(num,3:end-2,3:end-2).*reshape(duz_dxz,[1,nz,nx]) +co55(num,3:end-2,3:end-2).* reshape((dux_dzz + duz_dxz),[1,nz,nx])+...
%         co11_dx(num,:,:).* reshape(dux_dx,[1,nz,nx])+ co13_dx(num,:,:) .* reshape(duz_dz,[1,nz,nx]) + co55_dz(num,:,:) .* reshape((dux_dz + duz_dx),[1,nz,nx]))+...
%         coe2(num,3:end-2,3:end-2).*(co11(num,3:end-2,3:end-2).*reshape(dux_dxxo,[1,nz,nx]) +...
%         co13(num,3:end-2,3:end-2).*reshape(duz_dxzo,[1,nz,nx]) +co55(num,3:end-2,3:end-2).* reshape((dux_dzzo + duz_dxzo),[1,nz,nx])+co11_dx(num,:,:).* reshape(dux_dxo,[1,nz,nx])+ co13_dx(num,:,:) .* reshape(duz_dzo,[1,nz,nx]) +...
%         co55_dz(num,:,:) .* reshape((dux_dzo + duz_dxo),[1,nz,nx]))+e(num).*ra(num,:,:); 
% 
%     rb(num,:,:)=coe1(num,3:end-2,3:end-2).*(co55(num,3:end-2,3:end-2).*reshape((dux_dxz + duz_dxx),[1,nz,nx]) + co13(num,3:end-2,3:end-2).*reshape(dux_dxz,[1,nz,nx]) + co33(num,3:end-2,3:end-2).*reshape(duz_dzz,[1,nz,nx]) +...
%         co55_dx(num,:,:).* reshape((dux_dz + duz_dx),[1,nz,nx])+ co33_dz(num,:,:).* reshape(duz_dz,[1,nz,nx]) + co13_dz(num,:,:) .* reshape(dux_dx,[1,nz,nx]))+...
%         coe2(num,3:end-2,3:end-2).*(co55(num,3:end-2,3:end-2).*reshape((dux_dxzo + duz_dxxo),[1,nz,nx]) +...
%         co13(num,3:end-2,3:end-2).*reshape(dux_dxzo,[1,nz,nx]) + co33(num,3:end-2,3:end-2).* reshape(duz_dzzo,[1,nz,nx])+ co55_dx(num,:,:).* reshape((dux_dzo + duz_dxo),[1,nz,nx])+...
%         co33_dz(num,:,:) .* reshape(duz_dzo,[1,nz,nx]) + co13_dz(num,:,:) .* reshape(dux_dxo,[1,nz,nx]))+e(num).*rb(num,:,:); 
%  
%  end 
%      ux3(3:end-2,3:end-2) = 2.0*ux2(3:end-2,3:end-2) - ux1(3:end-2,3:end-2) + (c11(3:end-2,3:end-2) .* dux_dxx + c13(3:end-2,3:end-2) .* duz_dxz +...
%         c55(3:end-2,3:end-2) .* (dux_dzz + duz_dxz) + c11_dx .* dux_dx + c13_dx .* duz_dz + c55_dz .* (dux_dz + duz_dx)+...
%         reshape(sum(ra),[nz,nx]) ).*dt2rho(3:end-2,3:end-2);
%     
%      uz3(3:end-2,3:end-2) = 2.0*uz2(3:end-2,3:end-2) - uz1(3:end-2,3:end-2)+(c55(3:end-2,3:end-2) .* (dux_dxz + duz_dxx) + c13(3:end-2,3:end-2) .* dux_dxz +...
%         c33(3:end-2,3:end-2) .* duz_dzz + c55_dx .* (dux_dz + duz_dx) +  c33_dz .* duz_dz + c13_dz .* dux_dx+...
%         reshape(sum(rb),[nz,nx]) ).*dt2rho(3:end-2,3:end-2);
%%    % Add source term
    ux3(jsrc, isrc) = ux3(jsrc, isrc) + force_x(it);
    uz3(jsrc, isrc) = uz3(jsrc, isrc) + force_z(it);
    % Exchange data between t-2 (1), t-1 (2) and t (3) and apply ABS
    ux1 = ux2 .* weights;
    ux2 = ux3 .* weights;
    uz1 = uz2 .* weights;
    uz2 = uz3 .* weights;
    % Output
    if mod(it,IT_DISPLAY) == 0
%         fprintf('Time step: %d \t %.4f s\n',it, single(t(it)));
%         u=sqrt(ux3.^2 + uz3.^2);
%         imagesc(u); colorbar; colormap jet; xlabel('m'); ylabel('m');
%         title(['Step = ',num2str(it),'/',num2str(nt),', Time: ',sprintf('%.4f',t(it)),' sec']);
%         axis equal tight; drawnow;
        snap_ux(:,:,i_snap)=ux3;
        snap_uz(:,:,i_snap)=uz3;
        i_snap=i_snap+1;
    end
    rec_x(it,1)=it*dt;
    rec_z(it,1)=it*dt;
    for rec_num=2:168
    rec_x(it,rec_num)=ux3(465,509+(rec_num-1)*10);
    rec_z(it,rec_num)=uz3(465,509+(rec_num-1)*10);
    end
    dux_dxxo=dux_dxx;
    dux_dxzo=dux_dxz;
    dux_dzzo=dux_dzz;
    duz_dxxo=duz_dxx;
    duz_dxzo=duz_dxz;
    duz_dzzo=duz_dzz;
    dux_dxo=dux_dx;
    dux_dzo=dux_dz;
    duz_dxo=duz_dx;
    duz_dzo=duz_dz;      
    ra_old=ra;
    rb_old=rb;
end
toc; disp('End');
        save('disc_x_m','snap_ux')
        save('disc_z_m','snap_uz')
        save('recec_x_m','rec_x')
        save('recec_z_m','rec_z')