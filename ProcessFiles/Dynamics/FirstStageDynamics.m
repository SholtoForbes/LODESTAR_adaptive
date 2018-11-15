function [dz,q,xi,vec_angle,T,D,dm] = FirstStageDynamics(z,u,t,phase,interp,Throttle,Stage1,Atmosphere,auxdata)
% function [dz,q,phi] = rocketDynamics(z,u,t,phase,scattered)
%This function determines the dynamics of the system as time derivatives,
%from input primals. Called by FirstStageDynamics


warning('off','MATLAB:griddedInterpolant:MeshgridEval2DWarnId'); % turn off warning ID which warns of slightly inefficient meshgrid usage

h = z(1,:);   %Height
v = z(2,:);   %Velocity
m = z(3,:);   %Mass
gamma = z(4,:);
alpha = z(5,:);
zeta = z(6,:);

dalphadt = z(7,:);

phi = z(8,:);

dalphadt2 = u(1,:); % control is second derivative of AoA with time

if isnan(gamma)
    gamma = 1.5708;
end

%%%% Compute gravity from inverse-square law:
rEarth = 6.3674447e6;  %(m) radius of earth
mEarth = 5.9721986e24;  %(kg) mass of earth
G = 6.67e-11; %(Nm^2/kg^2) gravitational constant
g = G*mEarth./((h+rEarth).^2);

g0 = 9.81;

if h>=0 & strcmp(phase,'postpitch')
density = interp1(Atmosphere(:,1),Atmosphere(:,4),h);
P_atm = interp1(Atmosphere(:,1),Atmosphere(:,3),h);
speedOfSound = interp1(Atmosphere(:,1),Atmosphere(:,5),h);
else
    density = interp1(Atmosphere(:,1),Atmosphere(:,4),0);
P_atm = interp1(Atmosphere(:,1),Atmosphere(:,3),0);
speedOfSound = interp1(Atmosphere(:,1),Atmosphere(:,5),0);
end

q = 0.5*density.*v.^2;
T = Stage1.TSL + (101325 - P_atm)*Stage1.TMod; % Thrust from Falcon 1 users guide. exit area calculated in SCALING.docx
T = T.*Throttle; % Throttle down

Isp = Stage1.IspSL + (101325 - P_atm)*Stage1.IspMod; % linear regression of SL and vacuum Isp. From encyclopaedia astronautica, backed up by falcon 1 users guide
dm = -T./Isp./g0;

mach = v./speedOfSound;

% interpolate coefficients
Cd = interp.DragGridded(mach,rad2deg(alpha));
Cl = interp.LiftGridded(mach,rad2deg(alpha));
Cm = interp.MomentGridded(mach,rad2deg(alpha));


% Compute Thrust Vector
vec_angle = -asin(Cm./(T/(auxdata.Stage1.CG - auxdata.Stage1.L)./q));

%%%% Compute the drag and lift:
Area = Stage1.Area; 

D = 0.5*Cd.*Area.*density.*v.^2*auxdata.Cdmod;

L = 0.5*Cl.*Area.*density.*v.^2;

%modify lift and thrust for vectoring
L = L + T.*sin(vec_angle);
T = T.*cos(vec_angle);

switch phase
    case 'prepitch'
    gamma = 1.5708*ones(1,length(h)); % Control Trajectory Angle 
    zeta = 0;
    case 'postpitch'
    %Do nothing
end

xi = 0; 

% This determines the dynamics of the system.
% Set up like this because xi is a quasi-forward sim instead of a primal




[dr,dxi,dphi,dgamma,dv,dzeta] = RotCoordsFirst(h+rEarth,xi,phi,gamma,v,zeta,L,D,T,m,alpha,phase);



switch phase
    case 'prepitch'
    dgamma = 0; % Control Trajectory Angle 
    dzeta = 0;
    case 'postpitch'
    %Do nothing
end

dz = [dr;dv;dm;dgamma;dalphadt;dzeta;dalphadt2;dphi;dxi];

if any(isnan(dz))
    disp('NaN Values Detected')
end

end