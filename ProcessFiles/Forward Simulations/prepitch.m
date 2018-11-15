function cost = prepitch(controls,m_f,scattered,Throttle,Vehicle,Atmosphere,auxdata)

% h0_prepitch = 0;  %Rocket starts on the ground
h0_prepitch = controls(1);
t_f = controls(2);
v0_prepitch = 0;  %Rocket starts stationary

gamma0_prepitch = deg2rad(90);

Isp = Vehicle.IspSL;
T = Vehicle.TSL*Throttle;
dm = -T./Isp./9.81;
m0_prepitch = m_f - dm*t_f;  %Add some fuel onto the rocket to go first section

phase = 'prepitch';
tspan = [0 t_f]; % time to fly before pitchover (ie. straight up)

y0 = [h0_prepitch, v0_prepitch, m0_prepitch, gamma0_prepitch, 0, 0, 0,0,0];

% this performs a forward simulation before pitchover. The end results of
% this are used as initial conditions for the optimiser. 
[t_prepitch, y] = ode45(@(t,y) FirstStageDynamics(y,0,0,phase,scattered,Throttle,Vehicle,Atmosphere,auxdata), tspan, y0);  

cost = (90 - y(end,1))^2 + (15 - y(end,2))^2;
