% Max Louyot, April 2019
clc, close all, clear all


% This code computes the frequencies (rad/s or Hz), modeshapes and NAVMI
% factor of a free-clamped ring rotating in a fluid for a given number of
% nodal diameters and circles; the geometry is customizable. It also allows
% for simple Dirac forced response analysis.
% This code refers to equations from Leissa (1969), Amabili et al. (1996)
%% Parameters

%~Output unit (Hz or rad/s)
hertz = true; % 'false' for rad/s, 'true' for Hz
% always use rad/s values for input!

%~Material
    E = 200*10^9; % Young's modulus (Pa)
    nu = 0.27; % Poisson's ratio
    rhoD = 7680; % mass density of the solid (kg/m3)
%~Fluid
    rhoF = 997; % mass density of the fluid (kg/m3)
    K = 0.4; % entrainment coefficient (see Poncet et al., 2005)

%~Geometry
    a = 0.2; % outer radius of the plate (m)
    b = 0.025; % inner radius of the plate (m)
    % always works for b/a <= 0.125, otherwise may need correction
    h = 0.008; % plate thickness (m)
    Hup = 0.01; % top axial gap (m)
    Hdown = 0.097; % bottom axial gap (m)

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

Ms = rhoD*h;                    % structural modal mass (kg)
Ks = Df*k1^4;                   % structural modal stiffness (N/m)

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
    +D*besselk(n,k1.*r)).*cos(n.*theta);    % displacement
figure(1);
surf(r.*cos(theta), r.*sin(theta), W);      % modeshape

%% Added mass M -- no stator coupling

syms eta rho    % rho = r/a

HA = @(eta) 1./(lambda^2-eta.^2).*(lambda.*besselj(n,eta).*besselj(n+1,...
    lambda)-eta.*besselj(n+1,eta).*besselj(n,lambda))-b/a./(lambda^2-...
    eta.^2).*(lambda.*besselj(n,b/a.*eta).*besselj(n+1,b/a*lambda)...
    -eta.*besselj(n+1,b/a.*eta).*besselj(n,b/a*lambda));
HB = @(eta) 1./(lambda^2-eta.^2).*(lambda.*besselj(n,eta).*bessely(n+1,...
    lambda)-eta.*besselj(n+1,eta).*bessely(n,lambda))-b/a./(lambda^2-...
    eta.^2).*(lambda.*besselj(n,b/a.*eta).*bessely(n+1,b/a*lambda)...
    -eta.*besselj(n+1,b/a.*eta).*bessely(n,b/a*lambda));
HC = @(eta) 1./(lambda^2+eta.^2).*(lambda.*besselj(n,eta).*besseli(n+1,...
    lambda)+eta.*besselj(n+1,eta).*besseli(n,lambda))-b/a./(lambda^2+...
    eta.^2).*(lambda.*besselj(n,b/a.*eta).*besseli(n+1,b/a*lambda)...
    +eta.*besselj(n+1,b/a.*eta).*besseli(n,b/a*lambda));
HD = @(eta) 1./(lambda^2+eta.^2).*(-lambda.*besselj(n,eta).*besselk(n+1,...
    lambda)+eta.*besselj(n+1,eta).*besselk(n,lambda))-b/a./(lambda^2+...
    eta.^2).*(-lambda.*besselj(n,b/a.*eta).*besselk(n+1,b/a*lambda)...
    +eta.*besselj(n+1,b/a.*eta).*besselk(n,b/a*lambda));
H = @(eta) A.*HA(eta)+B.*HB(eta)+C.*HC(eta)+D.*HD(eta);
    % obtained from resolution with Hankel transform

W = @(rho) (A.*besselj(n,lambda.*rho)+B.*bessely(n,lambda.*rho)...
    +C.*besseli(n,lambda.*rho)+D.*besselk(n,lambda.*rho));

fun = @(eta,rho) rhoF*a*H(eta).*besselj(n,eta.*rho).*W(rho);
M = integral2(fun,0,Inf,b/a,1,'AbsTol',1e-2)     % added mass (kg)
% expression is obtained through the calculation of the work rate and
% the integration on the top and bottom surfaces of the disk
% NB: the rotation does not change the value of M

%% NAVMI factor -- no stator coupling

syms eta rho    % rho = r/a

if n == 0       % integral of cos(n*theta)^2 over 0..2*pi
    psi = 2*pi;
else
    psi = pi;
end

Tf_d = @(eta,rho) a^3*psi*rhoF/2*rho.*W(rho).*H(eta).*besselj(n,...
    eta.*rho).*(1+exp(-Hdown*2/a.*eta))./(1-exp(-Hdown*2/a.*eta));
    % kinetic energy of the fluid under the disk
Tf_u = @(eta,rho) a^3*psi*rhoF/2*rho.*W(rho).*H(eta).*besselj(n,...
    eta.*rho).*(1+exp(-Hup*2/a.*eta))./(1-exp(-Hup*2/a.*eta));
    % kinetic energy of the fluid above the disk
Tf = @(eta,rho) Tf_d(eta,rho)+Tf_u(eta,rho);
    % total kinetic energy of the fluid
Td = @(rho) a^2*psi*rhoD*h/2*rho.*W(rho).^2;
    % kinetic energy of the disk

beta0 = integral2(Tf,0,Inf,b/a,1,'AbsTol',1e-2)/integral(Td,b/a,1,...
    'AbsTol',1e-4);  % Added Virtual Mass Incremental factor -- no rotation
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

%% OmegaD range plot
% plots analytical model results over a specified rotation speed range

% beta0 = 1;      % artificially redefines beta0 (for regime test)
nb = 51;        % number of ploted points
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
forcedSin = figure(6);
hold on
plot(Omega(1,:),U,'k','linewidth',2)
xlabel('Excitation frequency $\Omega$ (Hz)','interpreter','latex',...
    'fontsize',14)
ylabel('Displacement amplitude $U$ (m)','interpreter','latex',...
    'fontsize',14)
title('Forced response ($\xi=0.1$)','interpreter','latex','fontsize',14)
set(gca,'ticklabelinterpreter','latex','fontsize',14)

OmegaR = Omega(1,find(U == max(U))) % max. amplitude resonance freq. (Hz)

