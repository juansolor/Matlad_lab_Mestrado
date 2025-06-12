function [sys,x0,str,ts]=BcsPlant(t,x,u,flag)
switch flag
    case 0
        [sys,x0,str,ts]=mdlInitializeSizes;
    case 1
        sys=mdlDerivatives(t,x,u);
    case 3
        sys=mdlOutputs(t,x,u);
    case {2, 4, 9 }
        sys = [];
    otherwise
        error(['Unhandled flag = ',num2str(flag)]);
end
function [sys,x0,str,ts]=mdlInitializeSizes
sizes = simsizes;
sizes.NumContStates = 5;
sizes.NumDiscStates = 0;
sizes.NumOutputs = 7;
sizes.NumInputs = 4;
sizes.DirFeedthrough = 0;
sizes.NumSampleTimes = 0;
sys=simsizes(sizes);
% x0=[7264733.79612223;672871.507300900;0.0123778175929964;50;90]';
x0 = [7427884.28284910;1051972.80242121;0.0119993084637901;50;65]';
str=[];
ts=[];

function sys=mdlDerivatives(t,x,u)
g   = 9.81;   % Gravitational acceleration constant [m/s²]
Cc = 2e-5 ;   % Choke valve constant
A1 = 0.008107;% Cross-section area of pipe below ESP [m²]
A2 = 0.008107;% Cross-section area of pipe above ESP [m²]
D1 = 0.1016;  % Pipe diameter below ESP [m]
D2 = 0.1016;  % Pipe diameter above ESP [m]
h1 = 200;     % Heigth from reservoir to ESP [m]
hw = 1000;    % Total vertical distance in well [m]
L1 =  500;    % Length from reservoir to ESP [m]
L2 = 1200;    % Length from ESP to choke [m]
V1 = 4.054;   % Pipe volume below ESP [m3]
V2 = 9.729;   % Pipe volume above ESP [m3]
f0 = 60;      % ESP characteristics reference freq [Hz]
q0_dt = 0.00542847756833871; % Downtrhust flow at f0 m3/s
q0_ut = 0.0134581949664923;  % Uptrhust flow at f0 m3/s
Inp = 65;     % ESP motor nominal current [A]
Pnp = 1.625e5;% ESP motor nominal Power [W]
b1 = 1.5e9;   % Bulk modulus below ESP [Pa]
b2 = 1.5e9;   % Bulk modulus above ESP [Pa]
M  = 1.992e8; % Fluid inertia parameters [kg/m4]
rho = 950;    % Density of produced fluid [kg/mÃ?Â³]
PI = 2.32e-9; % Well productivy index [m3/s/Pa]
mu  = 0.025;  % Viscosity [Pa*s]
dfq_max = 0.5;   % máxima variação em f/s
dzc_max = 1;   % máxima variação em zc %/s
tp =[1/dfq_max;1/dzc_max];  % Actuator Response time
CH = -0.03*mu + 1;
Cq = 2.7944*mu^4 - 6.8104*mu^3 + 6.0032*mu^2 - 2.6266*mu + 1;
Cp = -4.4376*mu^4 + 11.091*mu^3 -9.9306*mu^2 + 3.9042*mu + 1;

pr = u(3);%1.26e7; % Reservoir pressure
pm = u(4);%2e5;%101325;%1atm fez o foco do sinal ficar mais centrado no envelope   nominal: 2e6;    % Manifold pressure
zcref= u(2);
fqref = u(1);

%states 
pbh = x(1);
pwh = x(2);
q = x(3);
fq = x(4);
zc = x(5);


% Calculo do HEAD e delta de pressão
q0 = q/Cq*(f0/fq);
H0 = -1.2454e6*q0^2 + 7.4959e3*q0 + 9.5970e2;
H = CH*H0*(fq/f0)^2; % Head
Dp = rho*g*H;       % Delta de pressão
% Calculo da Potencia e corrente da bomba
P0 = -2.3599e9*q0^3 -1.8082e7*q0^2 +4.3346e6*q0 + 9.4355e4;
P = Cp*P0*(fq/f0)^3; % Potencia
I = Inp*P/Pnp;       % Corrente
% Calculo da pressão de intaike
F1 = 0.158*((rho*L1*q^2)/(D1*A1^2))*(mu/(rho*D1*q))^(1/4);
F2 = 0.158*((rho*L2*q^2)/(D2*A2^2))*(mu/(rho*D2*q))^(1/4);
pin = pbh - rho*g*h1 - F1;
% Vazao do rezervatorio vazao da chocke
qr  = PI*(pr - pbh);
qc  = Cc*(zc/100)*sign((pwh - pm))*sqrt(abs(pwh - pm));


sys(1) = b1/V1*(qr - q);%dpbhdt 
sys(2) = b2/V2*(q - qc);%dpwhdt
sys(3) = 1/M*(pbh - pwh - rho*g*hw - F1 - F2 + Dp);%dqdt
sys(4) = (fqref - fq)/tp(1);%dfqdt
sys(5) = (zcref - zc)/tp(2);%dzcdt


% sys(1)=-(b1*(x(3) - PI*(pr - x(1))))/V1;
% sys(2)=(b2*(x(3) + (Cc*x(5)*abs(pm - x(2))^(1/2)*sign(pm - x(2)))/100))/V2;
% sys(3)= -(x(2) - x(1) + g*hw*rho - (CH*g*rho*x(4)^2*((74959*f0*x(3))/(10*Cq*x(4)) - (1245400*f0^2*x(3)^2)/(Cq^2*x(4)^2) + 9597/10))/f0^2 + (79*L1*rho*x(3)^2*(mu/(D1*rho*x(3)))^(1/4))/(500*A1^2*D1) + (79*L2*rho*x(3)^2*(mu/(D2*rho*x(3)))^(1/4))/(500*A2^2*D2))/M;
% sys(4)=(u(1) - x(4))/tp(1);
% sys(5)= -(x(5) - zcref)/tp(2);


if abs(sys(4))>dfq_max
   sys(4) = dfq_max*sign(sys(4));
end
if abs(sys(5))>dzc_max
   sys(5) = dzc_max*sign(sys(5));
end

function sys=mdlOutputs(t,x,u)
g   = 9.81;   % Gravitational acceleration constant [m/s²]
Cc = 2e-5 ;   % Choke valve constant
A1 = 0.008107;% Cross-section area of pipe below ESP [m²]
A2 = 0.008107;% Cross-section area of pipe above ESP [m²]
D1 = 0.1016;  % Pipe diameter below ESP [m]
D2 = 0.1016;  % Pipe diameter above ESP [m]
h1 = 200;     % Heigth from reservoir to ESP [m]
hw = 1000;    % Total vertical distance in well [m]
L1 =  500;    % Length from reservoir to ESP [m]
L2 = 1200;    % Length from ESP to choke [m]
V1 = 4.054;   % Pipe volume below ESP [m3]
V2 = 9.729;   % Pipe volume above ESP [m3]
f0 = 60;      % ESP characteristics reference freq [Hz]
q0_dt = 0.00542847756833871; % Downtrhust flow at f0 m3/s
q0_ut = 0.0134581949664923;  % Uptrhust flow at f0 m3/s
Inp = 65;     % ESP motor nominal current [A]
Pnp = 1.625e5;% ESP motor nominal Power [W]
b1 = 1.5e9;   % Bulk modulus below ESP [Pa]
b2 = 1.5e9;   % Bulk modulus above ESP [Pa]
M  = 1.992e8; % Fluid inertia parameters [kg/m4]
rho = 950;    % Density of produced fluid [kg/mÃ?Â³]
PI = 2.32e-9; % Well productivy index [m3/s/Pa]
mu  = 0.025;  % Viscosity [Pa*s]
dfq_max = 0.5;   % máxima variação em f/s
dzc_max = 1;   % máxima variação em zc %/s
tp1= 1/dfq_max; % Actuator Response time
tp2=1/dzc_max;  % Actuator Response time
CH = -0.03*mu + 1;
Cq = 2.7944*mu^4 - 6.8104*mu^3 + 6.0032*mu^2 - 2.6266*mu + 1;
Cp = -4.4376*mu^4 + 11.091*mu^3 -9.9306*mu^2 + 3.9042*mu + 1;

pr = u(3);%1.26e7; % Reservoir pressure
pm = u(4);%2e5;%101325;%1atm fez o foco do sinal ficar mais centrado no envelope   nominal: 2e6;    % Manifold pressure
zcref=u(2);



q0 = x(3)/Cq*(f0/x(4));
H0 = -1.2454e6*q0^2 + 7.4959e3*q0 + 9.5970e2;
H = CH*H0*(x(4)/f0)^2; % Head

F1 = 0.158*((rho*L1*x(3)^2)/(D1*A1^2))*(mu/(rho*D1*x(3)))^(1/4);
pin = x(1) - rho*g*h1 - F1;

sys(1)=x(1);
sys(2)=x(2);
sys(3)=x(3);
sys(4)=x(4);
sys(5)=x(5);

 sys(6)= pin;%x(1) - g*h1*rho - (79*L1*rho*x(3)^2*(mu/(D1*rho*x(3)))^(1/4))/(500*A1^2*D1);%Pin
sys(7)=H;






