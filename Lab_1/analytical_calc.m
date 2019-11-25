clear var
close all

N = 60;
I = 5;
h = 0.0001:0.0000005:0.0004;
ag = 0.0003175;
wi = 0.0254;
mu0 = 4*pi*10^-7;
mur = 2000;
Lz = 0.00648;
we = wi/8;
wec = wi/4;
wc = wi/4;
hc = wi/2;

%% Magnetic flux density
B = Lz*we*(mu0*N*I*sqrt(2))./(2.*h + 2*wi/mur);
Bag = (mu0*N*I*sqrt(2))/(2*ag + 2*wi/mur);

figure
plot(h*1000,B)
xlabel('mm')
ylabel('T')

%% Inductance
% Air gap only
Rm1 = ag/(mu0*Lz*wec);
Rm2 = ag/(mu0*Lz*we);
Rm = Rm1 + Rm2/2;
L1 = N^2/Rm;

% With leakage
Rm3 = 3*wc/(mu0*Lz*hc);
Rm = (1/(Rm1 + Rm2/2) + 1/(Rm3/2))^-1;
L2 = N^2/Rm;

% With fringing
beta = wc/(2*ag);
alpha = 4*(beta*atan(beta)-log(sqrt(1+beta^2)))/pi;
kf = beta - alpha/2;
Rm1 = ag/(mu0*Lz*(wec+2*kf*ag));
Rm2 = ag/(mu0*Lz*(we+2*kf*ag));
Rm = (1/(Rm1 + Rm2/2) + 1/(Rm3/2))^-1;
L3 = N^2/Rm;

%% Attraction force
F = Bag^2*Lz*(2*we+wec)/(2*mu0);