% Max Louyot, May 2019
clc, close all, clear all

% WARNING: work in progress, results may be wrong!

% This code computes the frequencies (rad/s or Hz), modeshapes and NAVMI
% factor of a free-clamped ring rotating in a fluid for a given number of
% nodal diameters and circles; the geometry is customizable. It also allows
% for simple Dirac forced response analysis.
% This code refers to equations from Leissa (1969), Amabili et al. (1996)
% and Presas et al. (2015,2016)
%% Parameters

%~Output unit (Hz or rad/s)
hertz = true; % 'false' for rad/s, 'true' for Hz
geom = "Presas"; % 'ETH','Presas','custom' for automatic parameters input
% always use rad/s values for input!

if geom == "custom"
%~Material
    E = 200*10^9; % Young's modulus (Pa)
    nu = 0.27; % Poisson's ratio
    rhoD = 7680; % mass density of the solid (kg/m3)
%~Fluid
    rhoF = 0.997; % mass density of the fluid (kg/m3)
    K = 0.45; % entrainment coefficient (see Poncet et al., 2005)

%~Geometry
    a = 0.2; % outer radius of the plate (m)
    b = 0.025; % inner radius of the plate (m)
    c = 0.207; % radius to outer wall (m)
    % always works for b/a <= 0.125, otherwise may need correction
    h = 0.008; % plate thickness (m)
    Hup = 0.01; % top axial gap (m)
    Hdown = 0.097; % bottom axial gap (m)

% Quick def. option
elseif geom == "ETH"
    E=210e9; nu=0.3; rhoD=7850; rhoF=997; K=0.45; a=0.1; b=0.015; h=0.002;
    Hup = 0.0018; Hdown = a;
%     load("RotorDataETH.mat")
elseif geom == "Presas"
    E=200e9; nu=0.27; rhoD=7680; rhoF=997; K=0.45; a=0.2; b=0.025; h=0.008;
    Hup = 0.01; Hdown = 0.097; c=0.207;
%     load("RotorDataPres.mat")
else
    disp("Incorrect 'geom' input, please refer to line 15.")
    return
end

Df = E*h^3/12/(1-nu^2); % flexural rigidity (N.m)

%~Rotation speed
OmegaD = 25.1327; % disk angular velocity (rad/s)
if hertz == true
    OmegaD = OmegaD/2/pi;   % (Hz)
end

%~Mode
n = 3; % diametrical mode
s = 0; % circular mode

syms k

%% k calculation
% k^4 = rhoD*h*omega^2/Df

%//////////////~ /!\ works for b/a <= 0.125
if (n==0 || n==1)
    initk = (3.75*s+1.6)/a; % corrected guess value
else
    initk = (3.35*s+1.15*n)/a; % corrected guess value
end
%//////////////~

A1 = besselj(n,k*b); % coefficients A1...D4 found with Maple symbolic
B1 = bessely(n,k*b); %   resolution of the boundary conditions, namely:
C1 = besseli(n,k*b); %   W(b)=dW(b)/dr=0 (clamped inside)
D1 = besselk(n,k*b); %   Vr(a)=Mr(a)=0 (free outside)
A2 = besselj(n+1,k*b)*k*b-n*besselj(n,k*b);
B2 = bessely(n+1,k*b)*k*b-n*bessely(n,k*b);
C2 = -besseli(n+1,k*b)*k*b-n*besseli(n,k*b);
D2 = besselk(n+1,k*b)*k*b-n*besselk(n,k*b);
A3 = k^3*a^3*besselj(n+1,k*a)-n^2*nu*besselj(n,k*a)...
    -n^2*nu*besselj(n+1,k*a)*k*a+n^2*besselj(n,k*a)-n^3*besselj(n,k*a)...
    +k*a*besselj(n+1,k*a)*n^2-n*k^2*a^2*besselj(n,k*a)...
    +n^3*nu*besselj(n,k*a);
B3 = k^3*a^3*bessely(n+1,k*a)-n^2*nu*bessely(n,k*a)...
    -n^2*nu*bessely(n+1,k*a)*k*a+n^2*bessely(n,k*a)-n^3*bessely(n,k*a)...
    +k*a*bessely(n+1,k*a)*n^2-n*k^2*a^2*bessely(n,k*a)...
    +n^3*nu*bessely(n,k*a);
C3 = k^3*a^3*besseli(n+1,k*a)-n^2*nu*besseli(n,k*a)...
    +n^2*nu*besseli(n+1,k*a)*k*a+n^2*besseli(n,k*a)-n^3*besseli(n,k*a)...
    -k*a*besseli(n+1,k*a)*n^2+n*k^2*a^2*besseli(n,k*a)...
    +n^3*nu*besseli(n,k*a);
D3 = -k^3*a^3*besselk(n+1,k*a)-n^2*nu*besselk(n,k*a)...
    -n^2*nu*besselk(n+1,k*a)*k*a+n^2*besselk(n,k*a)-n^3*besselk(n,k*a)...
    +k*a*besselk(n+1,k*a)*n^2+n*k^2*a^2*besselk(n,k*a)...
    +n^3*nu*besselk(n,k*a);
A4 = -nu*besselj(n+1,k*a)*k*a-k^2*a^2*besselj(n,k*a)...
    +besselj(n+1,k*a)*k*a+nu*n*besselj(n,k*a)-n^2*nu*besselj(n,k*a)...
    +n^2*besselj(n,k*a)-n*besselj(n,k*a);
B4 = -nu*bessely(n+1,k*a)*k*a-k^2*a^2*bessely(n,k*a)...
    +bessely(n+1,k*a)*k*a+nu*n*bessely(n,k*a)-n^2*nu*bessely(n,k*a)...
    +n^2*bessely(n,k*a)-n*bessely(n,k*a);
C4 = nu*besseli(n+1,k*a)*k*a+k^2*a^2*besseli(n,k*a)...
    -besseli(n+1,k*a)*k*a+nu*n*besseli(n,k*a)-n^2*nu*besseli(n,k*a)...
    +n^2*besseli(n,k*a)-n*besseli(n,k*a);
D4 = -nu*besselk(n+1,k*a)*k*a+k^2*a^2*besselk(n,k*a)...
    +besselk(n+1,k*a)*k*a+nu*n*besselk(n,k*a)-n^2*nu*besselk(n,k*a)...
    +n^2*besselk(n,k*a)-n*besselk(n,k*a);
M = [A1 B1 C1 D1   % M . {An Bn Cn Dn}' = 0
     A2 B2 C2 D2   % k verifies det(M)=0
     A3 B3 C3 D3
     A4 B4 C4 D4];

k1 = vpasolve(det(M)==0,k,initk);
k1 = double(subs(k1))

lambda = k1*a;      % frequency parameter

%% Omega calculation

Ms = rhoD*h*pi*(a^2-b^2);       % structural modal mass (kg)
Ks = Df*k1^4*pi*(a^2-b^2);      % structural modal stiffness (N/m)

omegaA = sqrt(Ks/Ms)            % central frequency in vacuum (rad/s)
if hertz == true
    omegaA = omegaA/2/pi % (Hz)
end

%% Modeshape

syms BB CC DD % W = A*Jn+B*Yn+C*In+D*Kn
A = 1; % arbitrary choice to close the system of equations
[B,C,D] = vpasolve(subs(M(1:3,:),k,k1)*[1 BB CC DD]'==[0 0 0]',[BB,CC,DD]);
B = double(subs(B));    % M . {A B C D}' = 0
C = double(subs(C));
D = double(subs(D));

%% Displacement W

[r,theta] = meshgrid(b:(a-b)/50:a,0:pi/24:2*pi);
W = (A*besselj(n,k1.*r)+B*bessely(n,k1.*r)+C*besseli(n,k1.*r)...
    +D*besselk(n,k1.*r)).*cos(n.*theta); % displacement
figure(1);
surf(r.*cos(theta), r.*sin(theta), W); % modeshape

if n==3 && s==0 && geom=="ETH"   % FEM modeshapes comparison for mode (3,0)
    surf(r.*cos(theta+pi/6), r.*sin(theta+pi/6), W);
    hold on
    M1 = csvread('ModeshapesFEM/rotorZ_r008.csv',1,0);
    M2 = csvread('ModeshapesFEM/rotorZ_theta0.csv',1,0);
    plot3(M1(:,1),M1(:,2),M1(:,4).*562.5,'g','linewidth',2)
    plot3(M2(:,1),M2(:,2),M2(:,4).*562.5,'r','linewidth',2)
end

%% Flow potential

CC1 = ((-k1*c*besselk(n+1,k1*c)+n*besselk(n,k1*c))*(C*besseli(n+1,k1*a)*a...
    *k1-D*besselk(n+1,k1*a)*a*k1-A*besselj(n+1,k1*a)*k1*a-B*bessely(n+1,...
    k1*a)*a*k1+n*(A*besselj(n,k1*a)+C*besseli(n,k1*a)+B*bessely(n,k1*a)+D...
    *besselk(n,k1*a))))/((k1*a*(-k1*c*besselk(n+1,k1*c)+n*besselk(n,k1...
    *c))*besseli(n+1,k1*a)+k1*a*(besseli(n+1,k1*c)*k1*c+n*besseli(n,k1*c))...
    *besselk(n+1,k1*a)+(-besseli(n,k1*a)*c*k1*besselk(n+1,k1*c)-c...
    *besselk(n,k1*a)*k1*besseli(n+1,k1*c)+n*(besseli(n,k1*a)*besselk(n,k1...
    *c)-besseli(n,k1*c)*besselk(n,k1*a)))*n));

CC2 = -((C*besseli(n+1,k1*a)*a*k1-D*besselk(n+1,k1*a)*a*k1-A*besselj(n+1,...
    k1*a)*k1*a-B*bessely(n+1,k1*a)*a*k1+n*(A*besselj(n,k1*a)+C*besseli(n,...
    k1*a)+B*bessely(n,k1*a)+D*besselk(n,k1*a)))*(besseli(n+1,k1*c)*k1*c+n...
    *besseli(n,k1*c)))/((k1*a*(-k1*c*besselk(n+1,k1*c)+n*besselk(n,k1*c))...
    *besseli(n+1,k1*a)+k1*c*(k1*a*besselk(n+1,k1*a)-n*besselk(n,k1*a))...
    *besseli(n+1,k1*c)+n*(besseli(n,k1*c)*a*k1*besselk(n+1,k1*a)-besseli(...
    n,k1*a)*c*k1*besselk(n+1,k1*c)+n*(besseli(n,k1*a)*besselk(n,k1*c)...
    -besseli(n,k1*c)*besselk(n,k1*a)))));

%% NAVMI factor -- no stator coupling

syms r

Psi = @(r) A*besselj(n,k1.*r)+B*bessely(n,k1.*r)+C*besseli(n,k1.*r)...
    +D*besselk(n,k1.*r);

if n == 0       % integral of cos(n*theta)^2 over 0..2*pi
    psi = 2*pi;
else
    psi = pi;
end

Tf_1 = @(r) 0.5*psi*rhoF/k1*cot(k1*Hup).*Psi(r).^2.*r;
        % kinetic energy of the fluid above the disk
Tf_2 = @(r) 0.5*psi*rhoF/k1*cot(-k1*Hdown).*Psi(r).^2.*r;
        % kinetic energy of the fluid under the disk
Tf_3 = @(r) 0.5*psi*rhoF/k1*cot(k1*Hup)*(CC1*besseli(n,k1.*r)+CC2*...
    besselk(n,k1.*r)).*Psi(r).*r;
        % kinetic energy of the fluid above, in the radial gap
Tf_4 = @(r) 0.5*psi*rhoF/k1*cot(-k1*Hdown)*(CC1*besseli(n,k1.*r)+CC2*...
    besselk(n,k1.*r)).*Psi(r).*r;
        % kinetic energy of the fluid under, in the radial gap
Tf_disk = @(r) Tf_1(r)+Tf_2(r);
        % total kinetic energy of the fluid above and below the disk
Tf_gap = @(r) Tf_3(r)+Tf_4(r);
        % total kinetic energy of the fluid in the radial gap
Td = @(r) 0.5*psi*rhoD*h.*Psi(r).^2.*r;
        % kinetic energy of the disk

beta0 = (integral(Tf_disk,b,a)+integral(Tf_gap,a,c))/integral(Td,b,a);
        % Added Virtual Mass Incremental factor -- no rotation
Gamma0 = beta0*rhoD/rhoF*h/a;
        % Nondimensionalized Added Virtual Mass Incremental factor

omegaB = [0;0];
for i = 1:2     % i=1 : n, i=2 : -n
    omegaB(i) = (sqrt(omegaA^2*(beta0+1)-beta0*(n*(1-K)*OmegaD)^2)...
        -beta0*n*(1-K)*OmegaD)/(beta0+1);
            % solving omega/sqrt(1+AVMI) when AVMI depends on omega;
            % here beta0 is actually AVMI factor when OmegaD = 0.
    n = -n;
end
omegaB      % co-rotating and counter-rotating waves frequency

%% Fluid modes
% fluid modes in a (r,z) plan, for a given theta
% depending on the theta position, fluid is dragged to or pulled from the
% rotating disk (because of orthogonal mode and rotation)

syms r z
theta = 0;    % angle (rad)

% % this runs slow... do it once and save data, then just load data.
% phi1 = @(r,z) cos(n*theta)/k1*Psi(r)*(cot(k1*Hup)*cos(k1*z)+sin(k1*z));
%         % potential above disk
% phi2 = @(r,z) cos(n*theta)/k1*Psi(r)*(cot(-k1*Hdown)*cos(k1*z)+sin(k1*z));
%         % potential below disk
% phi3 = @(r,z) cos(n*theta)/k1*(CC1*besseli(n,k1*r)+CC2*besselk(n,k1*r))...
%     *(cot(k1*Hup)*cos(k1*z)+sin(k1*z));
%         % potential above radial gap
% phi4 = @(r,z) cos(n*theta)/k1*(CC1*besseli(n,k1*r)+CC2*besselk(n,k1*r))...
%     *(cot(-k1*Hdown)*cos(k1*z)+sin(k1*z));
%         % potential above radial gap
% 
% u1 = @(r,z) diff(phi1(r,z),r);  % radial velocity
% u2 = @(r,z) diff(phi2(r,z),r);
% u3 = @(r,z) diff(phi3(r,z),r);
% u4 = @(r,z) diff(phi4(r,z),r);
% 
% w1 = @(r,z) diff(phi1(r,z),z);  % vertical velocity
% w2 = @(r,z) diff(phi2(r,z),z);
% w3 = @(r,z) diff(phi3(r,z),z);
% w4 = @(r,z) diff(phi4(r,z),z);
% 
% U = zeros(108,183); W = zeros(108,183);
% for i = 1:98
%     for j = 1:176
%         U(i,j) = eval(subs(u2(r,z),[r z],[(j-1)/1000+b (i-1)/1000-Hdown]));
%         W(i,j) = eval(subs(w2(r,z),[r z],[(j-1)/1000+b (i-1)/1000-Hdown]));
%     end
%     for j = 177:183
%         U(i,j) = eval(subs(u4(r,z),[r z],[(j-1)/1000+b (i-1)/1000-Hdown]));
%         W(i,j) = eval(subs(w4(r,z),[r z],[(j-1)/1000+b (i-1)/1000-Hdown]));
%     end
% end
% for i = 99:108
%     for j = 1:176
%         U(i,j) = eval(subs(u1(r,z),[r z],[(j-1)/1000+b (i-1)/1000-Hdown]));
%         W(i,j) = eval(subs(w1(r,z),[r z],[(j-1)/1000+b (i-1)/1000-Hdown]));
%     end
%     for j = 177:183
%         U(i,j) = eval(subs(u3(r,z),[r z],[(j-1)/1000+b (i-1)/1000-Hdown]));
%         W(i,j) = eval(subs(w3(r,z),[r z],[(j-1)/1000+b (i-1)/1000-Hdown]));
%     end
% end
% 
% save VelocityField.mat U W
load VelocityField.mat    % (U,W) velocity field computed with code above

[r,z] = meshgrid(b:0.001:c,-Hdown:0.001:Hup);      % r,z coordinates

stream = figure(8);
quiver(r,z,U,W,'showarrowhead','on')     % plot speed vectors
hold on
plot([b a],[0 0],'-k','linewidth',2)
streamline(r,z,U,W,b,0.001)         % plot streamline
% streamline(r,z,U,W,a/2,-Hdown/2)
% streamline(r,z,U,W,4*a/5,Hup/2)
axis([b c -Hdown Hup])
xlabel('$r$ (m)','interpreter','latex','fontsize',14)
ylabel('$z$ (m)','interpreter','latex','fontsize',14)
set(gca,'ticklabelinterpreter','latex','fontsize',14)


%% OmegaD range plot
% plots analytical model results over a specified rotation speed range

beta0 = 1;      % artificially redefines beta0 (for regime test)
nb = 201;       % number of ploted points
step = 1;       % rotation speed steps (in Hz)

speed = zeros(1,nb);        % rotation speeds
omega_an = zeros(2,nb);     % analytical results
omega_ce = zeros(1,nb);     % central value

for i = 1:nb
    speed(1,i) = (i-1)*step;    % from 0 to (nb-1)*step vector
    for j = 1:2
        omega_an(j,i) = (sqrt(omegaA^2*(beta0+1)-beta0*(n*(1-K)*...
            speed(1,i))^2)-beta0*n*(-1)^(j-1)*(1-K)*speed(1,i))/(beta0+1);
    end
    omega_ce(1,i) = mean(omega_an(:,i));
end

rangeO = figure(4);
hold on
Aplot = plot(speed,omega_an,'-','color',[1 0.5 0],'linewidth',2);
Cplot = plot(speed,omega_ce,'-k','linewidth',2);
plot([0 (nb-1)*step],[omega_ce(1,1) omega_ce(1,1)],'--k')
xlabel('Rotor speed $\Omega_D$ (Hz)','interpreter','latex','fontsize',14)
ylabel('Frequency $\omega$ (Hz)','interpreter','latex','fontsize',14)
legend([Aplot(1),Cplot(1)],{'$\pm n$, Analytical','Central value'},...
    'interpreter','latex','fontsize',14)
% title('Aluminium','interpreter','latex','fontsize',14)
set(gca,'ticklabelinterpreter','latex','fontsize',14)


%% Flutter instability
% plots analytical model results for omega > 0 and < 0 (n > 0), both Re
% and Im parts --> allows the detection of the critical speed for flutter
% instability.

nb = 1001;       % number of ploted points
step = 1;       % rotation speed steps (in Hz)

speed = zeros(1,nb);        % rotation speeds
omega_re = zeros(2,nb);     % analytical results, real part
omega_im = zeros(2,nb);     % analytical results, imaginary part

for i = 1:nb                % non-dimensionnal
    speed(1,i) = (i-1)*step/omegaA;    % from 0 to (nb-1)*step vector
    for j = 1:2
        omega_re(j,i) = real((-1)^j*(sqrt((beta0+1)-beta0*(n*...
            (1-K)*speed(1,i))^2)-beta0*n*(-1)^(j-1)*(1-K)*speed(1,i))/...
            (beta0+1));
        omega_im(j,i) = imag((-1)^j*(sqrt((beta0+1)-beta0*(n*...
            (1-K)*speed(1,i))^2)-beta0*n*(-1)^(j-1)*(1-K)*speed(1,i))/...
            (beta0+1));
    end
end

rangeF = figure(5);
subplot(2,1,1);
hold on
Fplot3 = plot(speed,omega_re,'k-','linewidth',1);
ylabel('$\Re(\omega)/\omega_v$','interpreter','latex','fontsize',14)
set(gca,'ticklabelinterpreter','latex','fontsize',14,'xtick',[])
subplot(2,1,2);
hold on
plot(speed,omega_im,'k-','linewidth',1);
ylabel('$\Im(\omega)/\omega_v$','interpreter','latex','fontsize',14)
xlabel('Rotor speed $\Omega_D/\omega_v$','interpreter','latex','fontsize',14)
set(gca,'ticklabelinterpreter','latex','fontsize',14)
% legend([Fplot2(1),Fplot3(1),Fplot4(1)],{'$n=2$','$n=3$','$n=4$'},'interpreter','latex','fontsize',14)


%% Stability map
% non-dimensional critical speed (flutter) plot for any mode

nb = 10001;          % number of ploted points
step = 0.01;       % beta0 steps

beta00 = zeros(1,nb);       % AVMI factors
omega_cr = zeros(1,nb);     % critical speeds (flutter) [-]

for i = 1:nb
    beta00(1,i) = (i-1)*step;    % from 0 to (nb-1)*step vector
    omega_cr(1,i) = sqrt(1+1/beta00(1,i));
end

stability = figure(6);
Splot = loglog(omega_cr,beta00,'-k','linewidth',2);
xlabel('$n\Omega_R/\omega_v$','interpreter','latex','fontsize',14)
ylabel('$\beta_0$','interpreter','latex','fontsize',14)
% legend([Aplot(1),Cplot(1)],{'$\pm n$, Analytical','Central value'},...
%     'interpreter','latex','fontsize',14)
% title('Aluminium','interpreter','latex','fontsize',14)
set(gca,'ticklabelinterpreter','latex','fontsize',14)


%% Forced response -- Dirac excitation -- no fluid

m = rhoD*h*pi*(a^2-b^2);    % mass of disk (kg)

%//////////////~
r0 = a;                                % force application radius (m)
OmegaD = 10;                           % disk rotation freq. (Hz)
Omega = meshgrid(0:1:2*omegaA);        % excitation frequency (Hz)
xi = 0.1;                              % damping coefficient
%//////////////~

F0 = r0*(A*besselj(n,k1*r0)+B*bessely(n,k1*r0)+C*besseli(n,k1*r0)+...
    D*besselk(n,k1*r0));    % force amplitude (N)

U = F0/8/pi/m*((pi^2*(omegaA^2-(Omega(1,:)+n*OmegaD).^2).^2+(xi*omegaA...
    .*(Omega(1,:)+n*OmegaD)).^2).^(-1/2)+(pi^2*(omegaA^2-(Omega(1,:)...
    -n*OmegaD).^2).^2+(xi*omegaA.*(Omega(1,:)-n*OmegaD)).^2).^(-1/2));
        % displacement amplitude (m)

titl_save = "forcedN"+n+"_omegD"+OmegaD+"_xi01";
forcedSin = figure(7);
hold on
plot(Omega(1,:),U,'k','linewidth',2)
xlabel('Excitation frequency $\Omega$ (Hz)','interpreter','latex','fontsize',14)
ylabel('Displacement amplitude $U$ (m)','interpreter','latex','fontsize',14)
title('Forced response ($\xi=0.1$)','interpreter','latex','fontsize',14)
set(gca,'ticklabelinterpreter','latex','fontsize',14)

OmegaR = Omega(1,find(U == max(U))) % max. amplitude resonance frequency (Hz)

% "save as" command: '-depsc' for vector format / '-dpng' for image
% print(forcedSin,titl_save,'-dpng')
