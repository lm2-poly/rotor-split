% Max Louyot, April 2019
clc, close all, clear all


% This code computes the frequencies (rad/s or Hz), modeshapes and NAVMI
% factor of a clamped circular plate in a fluid for a given number of
% nodal diameters and circles; the geometry is customizable. It also allows
% for simple Dirac forced response analysis.
% This code refers to equations from Leissa (1969), Amabili et al. (1996)
%% Parameters

%~Output unit (Hz or rad/s)
hertz = true; % false for rad/s, true for Hz
% always use rad/s values for input!

%~Material
    E = 210*10^9; % Young's modulus (Pa)
    nu = 0.3; % Poisson's ratio
    rhoD = 7850; % mass density (kg/m3)
%~Fluid
    rhoF = 997; % mass density of the fluid (kg/m3)
    K = 0.4; % rotor entrainment coefficient (see Poncet et al., 2005)

%~Geometry
    a = 0.125; % outer radius of the plate (m)
    h = 0.0015; % plate thickness (m)
    Hup = 0.05*a; % top axial gap (m)

Df = E*h^3/12/(1-nu^2); % flexural rigidity (N.m)

%~Rotation speed
OmegaD = 25.1327; % rotor angular velocity (rad/s)
if hertz == true
    OmegaD = OmegaD/2/pi;   % (Hz)
end
    
%~Mode
n = 3; % diametrical mode
s = 0; % circular mode

syms k

%% k calculation

initk = (1.4*n+3.3*s+3.2)/a; % corrected asymptotic value

A1 = besselj(n,k*a); % coefficients A1...D4 found with Maple symbolic
C1 = besseli(n,k*a); %   resolution of the boundary conditions, namely:
A2 = besselj(n+1,k*a)*k*a-n*besselj(n,k*a); %   W(a)=dW(a)/dr=0
C2 = -besseli(n+1,k*a)*k*a-n*besseli(n,k*a); % (clamped outside)
M = [A1 C1;A2 C2];

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

syms CC % W = A*Jn+C*In
A = 1; % arbitrary choice to close the system of equations
C = vpasolve(subs(M(1,:),k,k1)*[1 CC]'==0,CC);
C = double(subs(C));

%% Displacement W

[r,theta] = meshgrid(0:a/50:a,0:pi/24:2*pi);
W = (A*besselj(n,k1.*r)+C*besseli(n,k1.*r)).*cos(n.*theta); % displacement
figure(1);
surf(r.*cos(theta), r.*sin(theta), W); % modeshape

%% Added mass M -- no rotor coupling

syms eta rho    % rho = r/a

HA = @(eta) 1./(lambda^2-eta.^2).*(lambda.*besselj(n,eta).*besselj(n+1,...
    lambda)-eta.*besselj(n+1,eta).*besselj(n,lambda))-a./(lambda^2-...
    eta.^2).*(lambda.*besselj(n,a.*eta).*besselj(n+1,a*lambda)...
    -eta.*besselj(n+1,a.*eta).*besselj(n,a*lambda));
HC = @(eta) 1./(lambda^2+eta.^2).*(lambda.*besselj(n,eta).*besseli(n+1,...
    lambda)+eta.*besselj(n+1,eta).*besseli(n,lambda))-a./(lambda^2+...
    eta.^2).*(lambda.*besselj(n,a.*eta).*besseli(n+1,a*lambda)...
    +eta.*besselj(n+1,a.*eta).*besseli(n,a*lambda));
H = @(eta) A.*HA(eta)+C.*HC(eta);
    % obtained from resolution with Hankel transform

W = @(rho) (A.*besselj(n,lambda.*rho)+C.*besseli(n,lambda.*rho));

fun = @(eta,rho) H(eta).*besselj(n,eta.*rho).*W(rho);
M = rhoF*a*integral2(fun,0,Inf,0,1,'AbsTol',1e-5)/2     % added mass (kg)
% expression is obtained through the calculation of the work rate and
% the integration on the bottom surface of the disk

%% NAVMI factor -- no rotor coupling

syms eta rho    % rho = r/a

if n == 0       % integral of cos(n*theta)^2 over 0..2*pi
    psi = 2*pi;
else
    psi = pi;
end

Tf_d = @(eta,rho) a^3*psi*rhoF/2*rho.*W(rho).*H(eta).*besselj(n,...
    eta.*rho).*(1+exp(-Hup*2/a.*eta))./(1-exp(-Hup*2/a.*eta));
    % kinetic energy of the fluid under the disk
Td = @(rho) a^2*psi*rhoD*h/2*rho.*W(rho).^2;
    % kinetic energy of the disk
beta0 = integral2(Tf_d,0,Inf,0,1,'AbsTol',1e-5)...
    /integral(Td,0,1,'AbsTol',1e-5);      % AVMI = Tf_d/Td
    %  Added Virtual Mass Incremental factor -- no rotation
Gamma0 = beta0*rhoD/rhoF*h/a;
    % Nondimensionalized Added Virtual Mass Incremental factor

omegaB = [0;0];
for i = 1:2     % i=1 : n, i=2 : -n
    omegaB(i) = (sqrt(omegaA^2*(beta0+1)-beta0*(n*K*OmegaD)^2)...
        -beta0*n*K*OmegaD)/(beta0+1);
            % solving omega/sqrt(1+AVMI) when AVMI depends on omega;
            % here beta0 is actually AVMI factor when OmegaD = 0.
    n = -n;
end
omegaB      % co-rotating and counter-rotating waves frequency

%% Forced response -- Dirac excitation -- no fluid

m = rhoD*h*pi*a^2;    % mass of disk (kg)

%//////////////~
r0 = a;                                % force application radius (m)
OmegaD = 0;                            % disk rotation freq. (Hz)
Omega = meshgrid(0:1:2*omegaA);        % excitation frequency (Hz)
xi = 0.1;                              % damping coefficient
%//////////////~

F0 = r0*(A*besselj(n,k1*r0)+C*besseli(n,k1*r0));    % force amplitude (N)

U = F0/2/m*(((omegaA^2-(Omega(1,:)+n*OmegaD).^2).^2+4*(xi*omegaA...
    .*(Omega(1,:)+n*OmegaD)).^2).^(-1/2)+((omegaA^2-(Omega(1,:)...
    -n*OmegaD).^2).^2+4*(xi*omegaA.*(Omega(1,:)-n*OmegaD)).^2).^(-1/2));
        % displacement amplitude (m)

titl_save = "forcedN"+n+"_omegD"+OmegaD+"_xi01";
forcedSin = figure(5);
hold on
plot(Omega(1,:),U,'k','linewidth',2)
xlabel('Excitation frequency $\Omega$ (Hz)','interpreter','latex',...
    'fontsize',14)
ylabel('Displacement amplitude $U$ (m)','interpreter','latex',...
    'fontsize',14)
title('Forced response ($\xi=0.1$)','interpreter','latex','fontsize',14)
set(gca,'ticklabelinterpreter','latex','fontsize',14)

OmegaR = Omega(1,find(U == max(U))) % max. amplitude resonance freq. (Hz)

