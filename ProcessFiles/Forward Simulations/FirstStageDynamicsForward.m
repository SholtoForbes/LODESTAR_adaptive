function [dz] = FirstStageDynamicsForward(z,zeta,alpha,phase,interp,Throttle,Vehicle,Atmosphere,auxdata)
global mach

h = z(1,:);   %Height
v = z(2,:);   %Velocity
m = z(3,:);   %Mass
gamma = z(4,:);
phi = z(5,:);

% 
dalphadt = 0;

if isnan(gamma)
    gamma = 1.5708;
end

%%%% Compute gravity from inverse-square law:
rEarth = 6.3674447e6;  %(m) radius of earth
mEarth = 5.9721986e24;  %(kg) mass of earth
G = 6.67e-11; %(Nm^2/kg^2) gravitational constant
g = G*mEarth./((h+rEarth).^2);

density = interp1(Atmosphere(:,1),Atmosphere(:,4),h);
P_atm = interp1(Atmosphere(:,1),Atmosphere(:,3),h);
speedOfSound = interp1(Atmosphere(:,1),Atmosphere(:,5),h);

q = 0.5*density.*v.^2;

% SCALE = 1.;
% SCALE = 1; %this is engine exit area scale
% Merlin 1C engine 

T = auxdata.Stage1.TSL + (101325 - P_atm)*auxdata.Stage1.TMod; % Thrust from Falcon 1 users guide. exit area calculated in SCALING.docx

T = T.*Throttle; % Throttle down

Isp = auxdata.Stage1.IspSL + (101325 - P_atm)*auxdata.Stage1.IspMod; % linear regression of SL and vacuum Isp. From encyclopaedia astronautica, backed up by falcon 1 users guide


dm = -T./Isp./g;


mach = v./speedOfSound;

Cd = interp.DragGridded(mach,rad2deg(alpha));
Cl = interp.LiftGridded(mach,rad2deg(alpha));
Cm = interp.MomentGridded(mach,rad2deg(alpha));


% Compute Thrust Vector

vec_angle = -asin(Cm./(T/(auxdata.Stage1.CG - auxdata.Stage1.L)./q));
%%%% Compute the drag:

Area = auxdata.Stage1.Area; 
D = 0.5*Cd.*Area.*density.*v.^2*auxdata.Cdmod;
global L
L = 0.5*Cl.*Area.*density.*v.^2;

%modify lift and thrust for vectoring
L = L + T.*sin(vec_angle);
T = T.*cos(vec_angle);

switch phase
    case 'prepitch'
    gamma = 1.5708*ones(1,length(h)); % Control Trajectory Angle 
    case 'postpitch'
    %Do nothing
end

xi = 0; 

[dr,dxi,dphi,dgamma,dv,dzeta] = RotCoordsFirst(h+rEarth,xi,phi,gamma,v,zeta,L,D,T,m,alpha,phase);

% dzeta
if isnan(dgamma)
dgamma = 0;
end

dz = [dr;dv;dm;dgamma;dphi;dzeta];

if any(isnan(dz))
    disp('NaN Values Detected')
end
% dz = [dr;dv;dm;dgamma];
end