
function [c,ceq] = sizing(m0,r,cap,AR,m_payload,opt)

%scale inputs
scaleFactors = [1000 1 .01 10 ];
m0  = m0  * scaleFactors(1);
r   = r   * scaleFactors(2);
cap = cap * scaleFactors(3);
AR  = AR  * scaleFactors(4);

%Inputs
m=m0;

%Mission
range = 100000; %Range [m]
vc = 65; %Cruise speed [m/s]
t_hover = 4*60; %Time to hover [s]

%Constants
specific_energy = 200*3600; %Fully integrated pack level [W.s/kg]
specific_power = 3500; %Fully integrated motor [W/kg]
wire_density = 0.19; %Fully integrated wire density [kg/m]
vtip = 340*0.65; %Maximum tip speed for noise
rho = 1.225; %Air density
g = 9.81; %Gravity
CD0 = 0.05; %Aircraft drag coefficient
osEff=0.7; %Oswald efficiency
toc = 0.15; %Wing thickness ratio
m_other_frac = .2; %Mass fraction to account for landing gear and fuse
cl = 0.6; %Cruise lift coefficient
LoD = cl/(CD0+cl^2/(pi*osEff*AR)); %Lift to drag ratio
e_motor = 0.9*0.95; %motor+controller efficiency
e_prop = 0.8; %propeller efficiency
cutoff = 2*2*pi;
control_margin = 0.2; %Margin on top of motor out power for control

%Wing geometry
Sref = 2*m*g/(rho*vc^2*cl);
wingspan = sqrt(Sref*AR);
chord = wingspan/AR;
thk = chord*toc; %Wing thickness
load_factor = 3.8*1.5; %Max design load factor
cmocl = 0.15; %Design ratio of moment to lift coefficient

%Material properties
E=70e9; %Young's modulus [Pa]
G=30e9; %Shear modulus [Pa]
density=2800; %Material density [kg/m^3]
stress=275e6; %Stress allowable [Pa]
tau=stress/sqrt(3); %Shear stress allowable [Pa]
cripple = 10; %cap thickness to width for crippling

%Vehicle properties
m_fixed = 100; %Non-passenger payload (accounts for avioincs, communications, etc)
n_motors = 8; %Number of motors
thrust_margin = n_motors/(n_motors-2)*(1+control_margin); %Thrust margin including control margin

%Hover analysis
solidity = 0.15; %Assumed rotor solidity
rotor_cd0 = 0.012; %Assumed rotor airfoil drag
thrust = m*g/n_motors; %Hover thrust per motor
A_disk = pi*r^2; %Rotor disk area per rotor

P_indD = (thrust_margin*thrust)^1.5 / sqrt(2*rho*A_disk); %Design induced power
P_proD = rho*A_disk*vtip^3*solidity*rotor_cd0/8; %Profile power

vhover = vtip/sqrt(thrust_margin); %Hover tip speed
P_indH = thrust^1.5 / sqrt(2*rho*A_disk); %Hover induced power
P_proH = rho*A_disk*vhover^3*solidity*rotor_cd0/8; %Profile power

P_motorD = P_indD+P_proD; %Design power
P_motorH = P_indH+P_proH; %Hover power

%Rotor mass
%m_blade = c_rotor*t_rotor*r*rho; %blade mass
%force = m_blade*0.5*r*omega^2;
%stress = force/(c_rotor*t_rotor) + 6*moment*t_rotor/(c_rotor*t_rotor^3);
omega = 1.1*vtip/r;
moment = 0.5*thrust_margin*thrust*r; %blade root moment
c_rotor = 0.5*solidity*A_disk/r;
t_rotor = sqrt(6*moment/(c_rotor*(stress - 0.5*rho*(r*omega)^2)));
m_blade = c_rotor*t_rotor*r*density; %blade mass

%Cruise analysis
drag = m*g/LoD;
P_cruise = drag*vc/e_prop+m*vc^3/range; %drag plus min acceleration to get to vc
t_cruise = range/vc;

%Propulsion system masses
m_motor = max(P_motorD,P_cruise/n_motors)/specific_power;
m_wire = 2*1.5*wingspan*wire_density;
m_rotor = 2*m_blade;

%Wing analysis
%Loads
shear = 0.5*load_factor*m*g;
moment = 2*shear*wingspan/(3*pi);
relief = (m_motor+m_rotor)*wingspan*(1/2+1/4) + wingspan*m_wire*(5/96); %moment due to motors+rotors+wire
moment = moment - relief*load_factor;
            
%Spar
%cap = sqrt(moment/(stress*cripple*thk));
m_spar = 2*cripple*cap^2*density;
I_cap = 0.5*cripple*cap^2*thk^2;

%Skins
torsion = shear*cmocl*chord; %Torque at section
Ae = 0.5*thk*chord; %Enclosed area
skin = torsion/(2*tau*Ae);
m_skin = 2*chord*skin*density;

%Shear web
%I = web*thk^3/12
%Q = 0.5*web*thk^2
%tau = shear(1)*Q/(I*web);
web = 6*shear/(tau*thk);
m_web = thk*web*density;

%Total wing mass
m_wing = (m_spar+m_skin+m_web)*wingspan;
m_wing = sum(m_wing);

%Energy
energy = P_motorH*t_hover*n_motors/e_motor+P_cruise*t_cruise/e_motor;
cRate = max([(1+control_margin)*P_motorH*n_motors P_cruise])/energy;

%Mass
m_battery = energy/specific_energy;
m = m_payload+m_fixed+n_motors*(m_motor+m_rotor)+m_battery+m_wire+m_wing;
m_other = m_other_frac*m;
m = m_other+m;

%Frequency check
M_eq = 33/140*m_wing+m_motor*(0.5^3+1^3);
freq = 1.732*sqrt(E*I_cap/(M_eq*(wingspan/2)^3));

%Stress check
sy = moment/(cap^2*cripple*thk);
  
switch opt
    case 0
        c=m;
        c=c/scaleFactors(1);
    case 1
        c=[cutoff-freq;sy-stress];
        ceq=m-m0;
        
        c(1)=c(1)/cutoff;
        c(2)=c(2)/stress;
        ceq=ceq/scaleFactors(1);
    case 2
        c=[m_fixed m_motor*n_motors m_rotor*n_motors m_battery m_wire m_wing m_other m_payload];
        ceq=[wingspan chord thk cap web skin freq sy thrust/A_disk cRate];
end