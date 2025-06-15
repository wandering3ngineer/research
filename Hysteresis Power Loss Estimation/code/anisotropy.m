T50=50;        % Temperature of 50 [K]
Hc50=104.698;  % Coercive field at 50[K] in [Oe] 
Ms50=41.79;    % Magnetic saturation at 50 [K] in [emu/g]    
T300=300;      % Temperature of 300 [K]
Hc300=5.787;   % Coercive field at 300[K] in [Oe]    
Ms300=34.99;   % Magnetic saturation at 300[K] in [emu/g]

pd=5320;                      % Density [Kg/m3]
Hc=Hc50*(1e3/(4*pi()))        % Coercive field in [A/m]
Ms=Ms50*pd                    % Magnetic saturation in [A/m]
B=0.5;                        % Random easy axis constant
tm=136;                       % Measurement time [s]
t0=1e-11;                     % Attempt time [s]
kb=1.38064852e-23;            % Boltzmann constant [J/K]
d=13;                         % Particle diameter [nm]
V=(4/3)*pi()*((d*1e-9)/2)^3   % Particle volume [m3]
T=T50;                        % Temperature [K]

Keff=linspace(0,1e10,100);    % Anisotropy [J/m3]

LHS=Hc;
RHS=((2.*B.*Keff)./Ms).*(1-((kb.*T.*log(tm./t0))./(Keff.*V)).^B); 

clf;
figure(1); hold on;
plot(Keff,RHS,'r-');
plot(Keff,LHS*ones(size(Keff)),'b-');
hold off;

% Results: Keff(50)=1.86e9 [J/m3], Keff(300)=8.883e7;
% Compute the Keff(0) and T0 values
Keff50=1.86e9;
Keff300=8.883e7;

T5=5;
Keff5=5.9e5;
T280=280;
Keff280=1.1e4;
T0=(Keff280.*T5-Keff5.*T280)./(Keff280-Keff5)
Keff0=Keff280./(1-T280/T0)

Tb=257;
KeffB=Keff0*(1-Tb/T0)

tN=t0.*2.718^((KeffB.*V)/(kb.*Tb))
fN=1/tN
